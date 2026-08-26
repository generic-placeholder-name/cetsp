#pragma once

#if (defined(CETSP_HASH_SET_ABSEIL) + \
     defined(CETSP_HASH_SET_UNORDERED) + \
     defined(CETSP_HASH_SET_ORDERED)) != 1
#error "Exactly one CETSP hash-set backend must be selected"
#endif

#if defined(CETSP_HASH_SET_ABSEIL)
#include <absl/container/flat_hash_set.h>
#elif defined(CETSP_HASH_SET_UNORDERED)
#include <unordered_set>
#elif defined(CETSP_HASH_SET_ORDERED)
#include <set>
#endif

#if defined(CETSP_HASH_SET_ABSEIL)
template<typename T>
using HashSet = absl::flat_hash_set<T>;
#elif defined(CETSP_HASH_SET_UNORDERED)
template<typename T>
using HashSet = std::unordered_set<T>;
#elif defined(CETSP_HASH_SET_ORDERED)
template<typename T>
using HashSet = std::set<T>;
#endif
