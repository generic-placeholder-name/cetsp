#pragma once

#include <filesystem>
#include <string>
#include <string_view>
#include <vector>

namespace cetsp::postprocess::detail {

class TemporaryDirectory {
public:
    explicit TemporaryDirectory(std::string_view prefix);
    ~TemporaryDirectory();

    TemporaryDirectory(const TemporaryDirectory&) = delete;
    TemporaryDirectory& operator=(const TemporaryDirectory&) = delete;
    TemporaryDirectory(TemporaryDirectory&&) = delete;
    TemporaryDirectory& operator=(TemporaryDirectory&&) = delete;

    [[nodiscard]] const std::filesystem::path& path() const noexcept {
        return path_;
    }

private:
    std::filesystem::path path_;
};

int runProcess(
    const std::filesystem::path& executable,
    const std::vector<std::string>& arguments);

} // namespace cetsp::postprocess::detail
