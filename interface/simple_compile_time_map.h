#ifndef SIMPLE_COMPILE_TIME_MAP_H__
#define SIMPLE_COMPILE_TIME_MAP_H__

#include <utility> // pair
#include <array>

template<typename K, typename V>
using KeyValue = std::pair<K, V>;


template<typename K, typename V, std::size_t N>
constexpr auto mapSize(const std::array<KeyValue<K, V>, N>) {
  return N;
}

template<typename K, typename V, std::size_t N, typename CF>
static constexpr V getValue(const std::array<KeyValue<K, V>, N> map, const K key, size_t range, CF comp)
{
  return (range == 0) ? throw "Key not found" :
    comp(map[range - 1].first, key) ? map[range - 1].second :
    getValue(map, key, range - 1, comp);
}

/**
 * Get the first key for the desired value
 */
template<typename K, typename V, std::size_t N, typename CF>
static constexpr K getFirstKey(const std::array<KeyValue<K, V>, N> map, const V value, size_t range, CF comp)
{
  // going from top to bottom here, to make sure to pick up the first one
  return (range == N) ? throw "Value not found" :
    comp(map[range].second, value) ? map[range].first :
    getFirstKey(map, value, range + 1, comp);
}


using ParameterIndex = KeyValue<const char*, int>;

static constexpr bool StrComp(const char* a, const char* b)
{
  return (*a && *b) ? (*a == *b && StrComp(a + 1, b + 1)) : (!*a && !*b);
}

template<std::size_t N>
static constexpr int getParIdx(std::array<ParameterIndex, N> params, const char* name)
{
  return getValue(params, name, mapSize(params), StrComp);
}

static constexpr bool IntComp(const int a, const int b) { return a == b; }

template<std::size_t N>
static constexpr const char* getParName(std::array<ParameterIndex, N> params, const int index)
{
  return getFirstKey(params, index, 0, IntComp);
}


#endif

