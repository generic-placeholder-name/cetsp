#pragma once

#if (defined(CETSP_HASH_SET_ABSEIL) + \
     defined(CETSP_HASH_SET_UNORDERED) + \
     defined(CETSP_HASH_SET_ORDERED)) != 1
#error "Exactly one CETSP set backend must be selected"
#endif

#if defined(CETSP_HASH_SET_ABSEIL)
#include <absl/container/flat_hash_set.h>
#elif defined(CETSP_HASH_SET_UNORDERED)
#include <unordered_set>
#elif defined(CETSP_HASH_SET_ORDERED)
#include <set>
#endif

#if defined(CETSP_HASH_SET_ABSEIL)
template<typename Id>
using IdSet = absl::flat_hash_set<Id>;
#elif defined(CETSP_HASH_SET_UNORDERED)
template<typename Id>
using IdSet = std::unordered_set<Id>;
#elif defined(CETSP_HASH_SET_ORDERED)
template<typename Id>
using IdSet = std::set<Id>;
#endif
