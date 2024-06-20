#pragma once

#if defined(_OPENMP)

#if defined(__GNUC__)
template<typename T, typename U>
T fetch_and_add(T &x, U inc) {
    return __atomic_fetch_add(&x, inc, __ATOMIC_SEQ_CST);
}
template<typename T>
bool compare_and_swap(T &x, const T &old_val, const T &new_val) {
  return __sync_bool_compare_and_swap(&x, old_val, new_val);
}
#elif defined(_MSC_BUILD)
template<typename T, typename U>
T fetch_and_add(T &x, U inc) {
  return _InterlockedExchangeAdd(reinterpret_cast<long*>(&x), inc);
}
template<typename T>
bool compare_and_swap(T &x, const T &old_val, const T &new_val) {
  return _InterlockedCompareExchange(reinterpret_cast<long*>(&x), new_val, old_val) == old_val;
}
#else
#error "error: more platform atomic definitions needed"
#endif

#else
template<typename T, typename U>
T fetch_and_add(T &x, U inc) {
  T orig_val = x;
  x += inc;
  return orig_val;
}
template<typename T>
bool compare_and_swap(T &x, const T &old_val, const T &new_val) {
  if (x == old_val) {
    x = new_val;
    return true;
  }
  return false;
}
#endif

