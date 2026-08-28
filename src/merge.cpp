#include "merge.hpp"

#include "circle_geometry.hpp"
#include "debug.hpp"
#include "id_set.hpp"

#include <boost/container/small_vector.hpp>
#include <boost/geometry.hpp>
#include <boost/geometry/index/rtree.hpp>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iterator>
#include <limits>
#include <numbers>
#include <numeric>
#include <optional>
#include <random>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

namespace {

using BoxValue = std::pair<Box, TreeNodeId>;
using HeapEntry = std::pair<double, TreeNodeId>;

template<std::size_t Capacity>
using BoxValues = boost::container::small_vector<BoxValue, Capacity>;

void rotatePoint(Point& point, double cosine, double sine) {
    const double x = bg::get<0>(point);
    const double y = bg::get<1>(point);
    bg::set<0>(point, x * cosine - y * sine);
    bg::set<1>(point, x * sine + y * cosine);
}

} // namespace

void removeCoveringCircles(std::vector<Circle>& circles) {
    const std::size_t circleCount = circles.size();
    std::vector<std::size_t> order(circleCount);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
        return circles[a].r > circles[b].r;
    });

    bgi::rtree<BoxValue, bgi::rstar<16>> circleIndex;
    std::vector<bool> removed(circleCount, false);
    for (std::size_t id : order) {
        if (removed[id]) continue;
        BoxValues<16> candidates;
        circleIndex.query(
            bgi::covers(circles[id].center),
            std::back_inserter(candidates));
        for (const BoxValue& candidate : candidates) {
            Circle& covering = circles[candidate.second];
            if (bg::distance(covering.center, circles[id].center) +
                    circles[id].r <=
                covering.r) {
                removed[candidate.second] = true;
                circleIndex.remove(
                    std::make_pair(circleBox(covering), candidate.second));
            }
        }
        circleIndex.insert(std::make_pair(circleBox(circles[id]), id));
    }

    std::vector<Circle> survivors;
    survivors.reserve(circleCount);
    for (std::size_t id = 0; id < circles.size(); ++id) {
        if (!removed[id]) survivors.push_back(std::move(circles[id]));
    }
    circles = std::move(survivors);

    DBG("Done removing");
    DBG("Remaining circles: " << circles.size());
    for (std::size_t id = 0; id < circles.size(); ++id) {
        const Circle& circle = circles[id];
        DBG("Circle ID: " << id
            << ", Center: (" << bg::get<0>(circle.center) << ", "
            << bg::get<1>(circle.center) << ")"
            << ", Radius: " << circle.r);
    }
}

namespace {

class MergeTreeBuilder {
public:
    MergeTreeBuilder(
        const std::vector<Circle>& circles,
        std::mt19937_64& randomEngine);

    [[nodiscard]] std::vector<TreeNode> build();

private:
    struct Nearest {
        TreeNodeId id;
        double gap;
    };

    struct ClusterState {
        Circle circle;
        bool active = true;
        std::optional<Nearest> nearest;
        IdSet<TreeNodeId> dependents;

        explicit ClusterState(Circle value) : circle(std::move(value)) {}
    };

    struct MergeCandidate {
        TreeNodeId first;
        TreeNodeId second;
        double gap;
    };

    template<std::size_t CandidateCount>
    [[nodiscard]] std::optional<Nearest> findNearest(
        TreeNodeId treeNodeId) const;
    void setNearest(TreeNodeId treeNodeId, Nearest nearest);
    void clearNearest(TreeNodeId treeNodeId, bool eraseHeapEntry);
    [[nodiscard]] MergeCandidate popClosestPair();
    void mergeClosestPair(std::size_t mergeIndex);
    void assertValid() const;

    std::size_t leafCount_;
    std::mt19937_64& randomEngine_;
    double rotationCosine_;
    double rotationSine_;
    std::vector<ClusterState> clusters_;
    std::vector<TreeNode> treeNodes_;
    bgi::rtree<BoxValue, bgi::rstar<16>> circleIndex_;
    std::set<HeapEntry> mergeQueue_;
    std::size_t activeCount_ = 0;
};

MergeTreeBuilder::MergeTreeBuilder(
    const std::vector<Circle>& circles,
    std::mt19937_64& randomEngine)
    : leafCount_(circles.size()),
      randomEngine_(randomEngine),
      rotationCosine_(0.0),
      rotationSine_(0.0),
      activeCount_(circles.size()) {
    assert(circles.size() >= 2);
    clusters_.reserve(2 * leafCount_ - 1);
    treeNodes_.reserve(2 * leafCount_);

    std::uniform_real_distribution<> angleDistribution(
        0.0, 2 * std::numbers::pi);
    const double angle = angleDistribution(randomEngine_);
    rotationCosine_ = std::cos(angle);
    rotationSine_ = std::sin(angle);

    for (const Circle& input : circles) {
        Circle rotated = input;
        rotatePoint(
            rotated.center, rotationCosine_, rotationSine_);
        clusters_.emplace_back(rotated);
        treeNodes_.push_back(
            TreeNode::leaf(rotated.center, rotated.r));
    }

    for (TreeNodeId id = 0; id < leafCount_; ++id) {
        circleIndex_.insert({circleBox(clusters_[id].circle), id});
    }

    for (TreeNodeId id = 0; id < leafCount_; ++id) {
        const std::optional<Nearest> nearest = findNearest<16>(id);
        if (!nearest) {
            throw std::logic_error(
                "failed to find a nearest neighbor while building merge tree");
        }
        setNearest(id, *nearest);
    }

    DBG("Done initializing NNs");
    DBG("Nearest neighbors of each node:");
    for (TreeNodeId id = 0; id < leafCount_; ++id) {
        DBG("Node " << id
            << " -> NN: " << clusters_[id].nearest->id
            << " (gap: " << clusters_[id].nearest->gap << ")");
    }
    assertValid();
}

template<std::size_t CandidateCount>
std::optional<MergeTreeBuilder::Nearest> MergeTreeBuilder::findNearest(
    TreeNodeId treeNodeId) const {
    const Circle& circle = clusters_[treeNodeId].circle;
    BoxValues<CandidateCount> candidates;
    circleIndex_.query(
        bgi::nearest(circle.center, CandidateCount),
        std::back_inserter(candidates));

    std::optional<Nearest> nearest;
    double bestGap = std::numeric_limits<double>::infinity();
    for (const BoxValue& candidateValue : candidates) {
        const TreeNodeId candidateId = candidateValue.second;
        if (candidateId == treeNodeId) continue;
        const Circle& candidate = clusters_[candidateId].circle;
        const double gap = gapDist(
            circle.center,
            candidate.center,
            circle.r,
            candidate.r);
        if (gap < bestGap) {
            bestGap = gap;
            nearest = Nearest{candidateId, gap};
        }
    }
    return nearest;
}

void MergeTreeBuilder::setNearest(
    TreeNodeId treeNodeId,
    Nearest nearest) {
    ClusterState& cluster = clusters_[treeNodeId];
    if (!cluster.active || cluster.nearest) {
        throw std::logic_error(
            "nearest neighbor can only be set for an unqueued active cluster");
    }
    if (nearest.id >= clusters_.size() ||
        !clusters_[nearest.id].active ||
        nearest.id == treeNodeId) {
        throw std::logic_error(
            "nearest-neighbor link does not name another active cluster");
    }

    cluster.nearest = nearest;
    mergeQueue_.insert({nearest.gap, treeNodeId});
    clusters_[nearest.id].dependents.insert(treeNodeId);
}

void MergeTreeBuilder::clearNearest(
    TreeNodeId treeNodeId,
    bool eraseHeapEntry) {
    ClusterState& cluster = clusters_[treeNodeId];
    if (!cluster.nearest) {
        return;
    }

    const Nearest nearest = *cluster.nearest;
    if (eraseHeapEntry) {
        const std::size_t erased =
            mergeQueue_.erase({nearest.gap, treeNodeId});
        if (erased != 1) {
            throw std::logic_error(
                "active cluster's merge-queue entry could not be removed");
        }
    }
    clusters_[nearest.id].dependents.erase(treeNodeId);
    cluster.nearest.reset();
}

MergeTreeBuilder::MergeCandidate MergeTreeBuilder::popClosestPair() {
    while (!mergeQueue_.empty()) {
        const auto entry = mergeQueue_.begin();
        const auto [gap, firstId] = *entry;
        mergeQueue_.erase(entry);

        const ClusterState& first = clusters_[firstId];
        if (!first.active || !first.nearest ||
            first.nearest->gap != gap) {
            continue;
        }
        const TreeNodeId secondId = first.nearest->id;
        if (secondId >= clusters_.size() ||
            !clusters_[secondId].active) {
            throw std::logic_error(
                "selected nearest neighbor is not active");
        }
        return {firstId, secondId, gap};
    }

    throw std::logic_error(
        "merge queue became empty before the merge tree was complete");
}

void MergeTreeBuilder::mergeClosestPair(std::size_t mergeIndex) {
    const MergeCandidate candidate = popClosestPair();
    const TreeNodeId firstId = candidate.first;
    const TreeNodeId secondId = candidate.second;
    ClusterState& first = clusters_[firstId];
    ClusterState& second = clusters_[secondId];
    if (!first.nearest || first.nearest->id != secondId ||
        !second.nearest) {
        throw std::logic_error(
            "selected merge pair has incomplete nearest-neighbor state");
    }

    IdSet<TreeNodeId> affected = first.dependents;
    affected.insert(second.dependents.begin(), second.dependents.end());
    affected.erase(firstId);
    affected.erase(secondId);

    const TreeNodeId newId = clusters_.size();
    DBG("Merge #" << (mergeIndex + 1)
        << ": " << firstId << " + " << secondId
        << " (gap=" << candidate.gap << ") -> newId=" << newId);

    const Circle firstCircle = first.circle;
    const Circle secondCircle = second.circle;
    const auto [combinedCenter, combinedRadius] = makeCombinedCircle(
        firstCircle.center,
        firstCircle.r,
        secondCircle.center,
        secondCircle.r,
        randomEngine_);
    Circle combined{combinedCenter, combinedRadius};

    treeNodes_.push_back(TreeNode::branch(
        firstId,
        secondId,
        candidate.gap,
        combined.center,
        combined.r));

    // The first entry was already popped; the second remains queued.
    clearNearest(firstId, false);
    clearNearest(secondId, true);
    circleIndex_.remove({circleBox(firstCircle), firstId});
    circleIndex_.remove({circleBox(secondCircle), secondId});
    first.active = false;
    second.active = false;
    first.dependents.clear();
    second.dependents.clear();
    activeCount_ -= 2;
    DBG("Removed old circles");

    clusters_.emplace_back(combined);
    circleIndex_.insert({circleBox(combined), newId});
    ++activeCount_;

    const std::optional<Nearest> newNearest = findNearest<2>(newId);
    if (!newNearest && mergeIndex != leafCount_ - 2) {
        throw std::logic_error(
            "merged circle has no neighbor before the final merge");
    }
    if (newNearest) {
        setNearest(newId, *newNearest);
    }
    DBG("  New circle " << newId
        << " nn=" << (newNearest
               ? std::to_string(newNearest->id)
               : std::string("<none>"))
        << " gap=" << (newNearest
               ? newNearest->gap
               : std::numeric_limits<double>::infinity()));

    for (TreeNodeId affectedId : affected) {
        if (affectedId >= clusters_.size() ||
            !clusters_[affectedId].active ||
            !clusters_[affectedId].nearest) {
            throw std::logic_error(
                "reverse-neighbor link does not name an active cluster");
        }
        clearNearest(affectedId, true);
        const std::optional<Nearest> nearest = findNearest<2>(affectedId);
        if (!nearest) {
            throw std::logic_error(
                "failed to recompute an active cluster's nearest neighbor");
        }
        setNearest(affectedId, *nearest);
        DBG("  Update NN for cluster " << affectedId
            << " -> nn=" << nearest->id
            << " gap=" << nearest->gap);
    }

    assertValid();
}

void MergeTreeBuilder::assertValid() const {
#ifndef NDEBUG
    assert(circleIndex_.size() == activeCount_);
    std::size_t observedActive = 0;
    for (TreeNodeId id = 0; id < clusters_.size(); ++id) {
        const ClusterState& cluster = clusters_[id];
        if (!cluster.active) {
            assert(!cluster.nearest);
            assert(cluster.dependents.empty());
            continue;
        }

        ++observedActive;
        if (cluster.nearest) {
            assert(cluster.nearest->id < clusters_.size());
            assert(clusters_[cluster.nearest->id].active);
            assert(cluster.nearest->id != id);
            assert(mergeQueue_.contains({cluster.nearest->gap, id}));
            assert(clusters_[cluster.nearest->id].dependents.contains(id));
        } else {
            assert(activeCount_ == 1);
        }
        for (TreeNodeId dependent : cluster.dependents) {
            assert(dependent < clusters_.size());
            assert(clusters_[dependent].active);
            assert(clusters_[dependent].nearest);
            assert(clusters_[dependent].nearest->id == id);
        }
    }
    assert(observedActive == activeCount_);

    for (const HeapEntry& entry : mergeQueue_) {
        const auto [gap, id] = entry;
        assert(id < clusters_.size());
        assert(clusters_[id].active);
        assert(clusters_[id].nearest);
        assert(clusters_[id].nearest->gap == gap);
    }

    for (const BoxValue& indexed : circleIndex_) {
        assert(indexed.second < clusters_.size());
        assert(clusters_[indexed.second].active);
        assert(bg::equals(
            indexed.first,
            circleBox(clusters_[indexed.second].circle)));
    }
#endif
}

std::vector<TreeNode> MergeTreeBuilder::build() {
    for (std::size_t mergeIndex = 0;
         mergeIndex < leafCount_ - 1;
         ++mergeIndex) {
        mergeClosestPair(mergeIndex);
    }

    for (TreeNode& node : treeNodes_) {
        rotatePoint(node.center, rotationCosine_, -rotationSine_);
    }

    DBG("Combine complete. Total tree nodes=" << treeNodes_.size());
    return std::move(treeNodes_);
}

} // namespace

std::vector<TreeNode> buildMergeTree(
    const std::vector<Circle>& circles,
    std::mt19937_64& randomEngine) {
    if (circles.empty()) {
        return {};
    }
    if (circles.size() == 1) {
        return {TreeNode::leaf(circles.front().center, circles.front().r)};
    }

    MergeTreeBuilder builder(circles, randomEngine);
    return builder.build();
}
