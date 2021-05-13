#include "toml/toml.hpp"

#include <array>
#include <cstring>
#include <iostream>
#include <map>
#include <string>
#include <type_traits>
#include <vector>

// Some boiler-plate reducing getter wrappers for toml11
namespace toml_h {

namespace is_array_like_impl {
template <typename T> struct is_array_like : std::false_type {};
template <typename T, std::size_t N>
struct is_array_like<std::array<T, N>> : std::true_type {};
template <typename... Args>
struct is_array_like<std::vector<Args...>> : std::true_type {};
} // namespace is_array_like_impl

// type trait to utilize the implementation type traits as well as decay the
// type
template <typename T> struct is_array_like {
  static constexpr bool const value =
      is_array_like_impl::is_array_like<typename std::decay<T>::type>::value;
};

namespace is_vector_impl {
template <typename T> struct is_vector : std::false_type {};
template <typename... Args>
struct is_vector<std::vector<Args...>> : std::true_type {};
template <typename... Args>
struct is_vector<std::initializer_list<Args...>> : std::true_type {};
} // namespace is_vector_impl

// type trait to utilize the implementation type traits as well as decay the
// type
template <typename T> struct is_vector {
  static constexpr bool const value =
      is_vector_impl::is_vector<typename std::decay<T>::type>::value;
};

namespace is_array_impl {
template <typename T> struct is_array : std::false_type {
  static constexpr size_t const len = 0;
};
template <typename T, std::size_t N>
struct is_array<std::array<T, N>> : std::true_type {
  static constexpr size_t const len = N;
};
} // namespace is_array_impl

// type trait to utilize the implementation type traits as well as decay the
// type
template <typename T> struct is_array {
  static constexpr bool const value =
      is_array_impl::is_array<typename std::decay<T>::type>::value;
  static constexpr size_t const len =
      is_array_impl::is_array<typename std::decay<T>::type>::len;
};

template <typename T>
struct find_is_specialized
    : std::conditional<std::is_floating_point<T>::value, std::true_type,
                       std::false_type>::type {};
template <typename T, size_t N>
struct find_is_specialized<std::array<T, N>> : std::true_type {};
template <typename... Args>
struct find_is_specialized<std::vector<Args...>> : std::true_type {};
template <> struct find_is_specialized<bool> : std::true_type {};

inline toml::value parse_card(std::string card_name = "") {
  toml::value card_parse_result;
  try {
    card_parse_result = toml::parse(card_name);
  } catch (std::runtime_error const &e) {
    std::cout << "[ERROR]: Failed to parse toml card: " << card_name
              << ", with error: " << e.what() << std::endl;
    abort();
  }
  return card_parse_result;
}

// ********* Prototypes

// Allows implicit cast for integer to floating point for vector elements
template <typename T = toml::value, typename Key = std::string>
typename std::enable_if<
    is_vector<T>::value &&
        std::is_floating_point<typename T::value_type>::value,
    T>::type
find(toml::value const &v, Key const &k);

// covers forwarding for vectors of non floating types
template <typename T = toml::value, typename Key = std::string>
typename std::enable_if<
    is_vector<T>::value &&
        !std::is_floating_point<typename T::value_type>::value,
    T>::type
find(toml::value const &v, Key const &k);

// Checks std::array sizes for floating types
template <typename T = toml::value, typename Key = std::string>
typename std::enable_if<
    is_array<T>::value && std::is_floating_point<typename T::value_type>::value,
    T>::type
find(toml::value const &v, Key const &k);

// Checks std::array sizes and forwards for non floating types
template <typename T = toml::value, typename Key = std::string>
typename std::enable_if<
    is_array<T>::value &&
        !std::is_floating_point<typename T::value_type>::value,
    T>::type
find(toml::value const &v, Key const &k);

// Allows certain string values to be interpreted as bools
template <typename T = toml::value, typename Key = std::string>
typename std::enable_if<std::is_same<T, bool>::value, T>::type
find(toml::value const &v, Key const &k);

// Allows implicit cast for integer to floating point
template <typename T = toml::value, typename Key = std::string>
typename std::enable_if<std::is_floating_point<T>::value, T>::type
find(toml::value const &v, Key const &k);

// default non specialized forwarder for single keys
template <typename T = toml::value, typename Key = std::string>
typename std::enable_if<!find_is_specialized<T>::value, T>::type
find(toml::value const &v, Key const &k);

// forwarder for vectors of keys
template <typename T = toml::value, typename Key = std::string>
T find_rec(toml::value const &v,
           typename std::enable_if<is_vector<Key>::value, Key>::type const &ks);

// ********* implementations

// Allows implicit cast for integer to floating point for vector elements
template <typename T = toml::value, typename Key = std::string>
typename std::enable_if<
    is_vector<T>::value &&
        std::is_floating_point<typename T::value_type>::value,
    T>::type
find(toml::value const &v, Key const &k) {
  T rtn;
  for (auto const &i : v.at(k).as_array()) {
    if (i.is_integer()) {
      rtn.push_back(typename T::value_type(i.as_integer()));
    } else {
      rtn.push_back(typename T::value_type(i.as_floating()));
    }
  }
  return rtn;
}

// covers forwarding for vectors of non floating types
template <typename T = toml::value, typename Key = std::string>
typename std::enable_if<
    is_vector<T>::value &&
        !std::is_floating_point<typename T::value_type>::value,
    T>::type
find(toml::value const &v, Key const &k) {
  return toml::find<T>(v, k);
}

// Checks std::array sizes for floating types
template <typename T = toml::value, typename Key = std::string>
typename std::enable_if<
    is_array<T>::value && std::is_floating_point<typename T::value_type>::value,
    T>::type
find(toml::value const &v, Key const &k) {

  if (v.at(k).as_array().size() != is_array<T>::len) {
    std::cout << "[ERROR]: When parsing key: " << k
              << ", expected array of length " << is_array<T>::len
              << " but read one of length: " << v.at(k).as_array().size()
              << std::endl;
    abort();
  }

  T rtn;
  size_t ctr = 0;
  for (auto const &i : v.at(k).as_array()) {
    if (i.is_integer()) {
      rtn[ctr] = typename T::value_type(i.as_integer());
    } else {
      rtn[ctr] = typename T::value_type(i.as_floating());
    }
    ctr++;
  }
  return rtn;
}

// Checks std::array sizes and forwards for non floating types
template <typename T = toml::value, typename Key = std::string>
typename std::enable_if<
    is_array<T>::value &&
        !std::is_floating_point<typename T::value_type>::value,
    T>::type
find(toml::value const &v, Key const &k) {

  if (v.at(k).as_array().size() != is_array<T>::len) {
    std::cout << "[ERROR]: When parsing key: " << k
              << ", expected array of length " << is_array<T>::len
              << " but read one of length: " << v.at(k).as_array().size()
              << std::endl;
    abort();
  }
  return toml::find<T>(v, k);
}

// Allows certain string values to be interpreted as bools
template <typename T = toml::value, typename Key = std::string>
typename std::enable_if<std::is_same<T, bool>::value, T>::type
find(toml::value const &v, Key const &k) {

  if (v.at(k).is_string()) {
    std::string const &str = v.at(k).as_string();

    if (str == "on") {
      return true;
    }
    if (str == "off") {
      return false;
    }
    if (str == "true") {
      return true;
    }
    if (str == "false") {
      return false;
    }
  }

  return toml::find<T>(v, k);
}

// Allows implicit cast for integer to floating point
template <typename T = toml::value, typename Key = std::string>
typename std::enable_if<std::is_floating_point<T>::value, T>::type
find(toml::value const &v, Key const &k) {

  if (v.at(k).is_integer()) {
    return T(v.at(k).as_integer());
  }

  return toml::find<T>(v, k);
}

// default non specialized forwarder for single keys
template <typename T = toml::value, typename Key = std::string>
typename std::enable_if<!find_is_specialized<T>::value, T>::type
find(toml::value const &v, Key const &k) {
  return toml::find<T>(v, k);
}

// forwarder for vectors of keys
template <typename T = toml::value, typename Key = std::string>
T find_rec(toml::value const &v, std::vector<Key> const &ks) {
  if (!v.contains(ks.front())) {
    std::cout << "[ERROR]: Key " << ks.front() << " not found." << std::endl;
    (void)v.at(ks.front());
  }

  if (ks.size() > 1) {
    return toml_h::find_rec<T>(v.at(ks.front()),
                               std::vector<Key>(ks.begin() + 1, ks.end()));
  }
  return toml_h::find<T>(v, ks.front());
}

//****** Other helper methods

template <typename T, typename Key = std::string>
T find_or_rec(toml::value const &v, std::vector<Key> const &ks, T const &def) {

  if (!v.contains(ks.front())) {
    return def;
  }

  if (ks.size() > 1) {
    return toml_h::find_or_rec<T>(v.at(ks.front()),
                                  std::vector<Key>(ks.begin() + 1, ks.end()));
  }
  return toml_h::find<T>(v, ks.front());
}

template <typename T, typename Key = std::string>
T find_or(toml::value const &v, Key const &k, T const &def) {

  if (!v.contains(k)) {
    return def;
  }
  return toml_h::find<T>(v, k);
}

// These should come first in case people try to pass map initializer lists as
// map-like literals, the std::vector<std::string> form does weird stuff to them
template <typename T>
T ensure_valid_string_option(std::string const &val,
                             std::map<std::string, T> const &options) {

  for (auto const &opt : options) {
    if (val == opt.first) {
      return opt.second;
    }
  }

  std::cout << "[ERROR]: Read: " << val
            << ", but only the following are valid options: " << std::endl;
  for (auto const &opt : options) {
    std::cout << "\t" << opt.first << std::endl;
  }
  abort();
}

template <typename Key = std::string, typename T>
T ensure_valid_string_option(toml::value const &v, Key const &k,
                             std::map<std::string, T> const &options) {

  return ensure_valid_string_option(toml_h::find<std::string>(v, k), options);
}

inline int ensure_valid_string_option(std::string const &val,
                                      std::vector<std::string> const &options) {

  for (size_t i = 0; i < options.size(); ++i) {
    if (val == options[i]) {
      return i;
    }
  }

  std::cout << "[ERROR]: Read: " << val
            << ", but only the following are valid options: " << std::endl;
  for (size_t i = 0; i < options.size(); ++i) {
    std::cout << "\t" << options[i] << std::endl;
  }
  abort();
}

template <typename Key = std::string>
int ensure_valid_string_option(toml::value const &v, Key const &k,
                               std::vector<std::string> const &options) {

  return ensure_valid_string_option(toml::find<std::string>(v, k), options);
}

template <typename T, typename Key = std::string>
void set_if_present(T &target, toml::value const &v, Key const &k) {

  if (!v.contains(k)) {
    return;
  }
  target = toml_h::find<T>(v, k);
}

template <typename T, typename Key = std::string>
void set_if_present_rec(T &target, toml::value const &v,
                        std::vector<Key> const &ks) {

  if (!v.contains(ks.front())) {
    return;
  }

  if (ks.size() > 1) {
    set_if_present_rec(target, v.at(ks.front()),
                       std::vector<Key>(ks.begin() + 1, ks.end()));
    return;
  }

  target = toml_h::find<T>(v, ks.front());
}

template <typename T, typename Key = std::string, typename OptCollection>
void set_option_if_present(T &target, toml::value const &v, Key const &k,
                           OptCollection const &options) {

  if (!v.contains(k)) {
    return;
  }
  std::string opt_val = toml_h::find<std::string>(v, k);
  target = toml_h::ensure_valid_string_option(opt_val, options);
}

template <typename T, typename Key = std::string, typename OptCollection>
void set_option_if_present_rec(T &target, toml::value const &v,
                               std::vector<Key> const &ks,
                               OptCollection const &options) {

  if (!v.contains(ks.front())) {
    return;
  }
  if (ks.size() > 1) {
    set_option_if_present_rec(target, v.at(ks.front()),
                              std::vector<Key>(ks.begin() + 1, ks.end()),
                              options);
    return;
  }
  std::string opt_val = toml_h::find<std::string>(v, ks.front());
  target = toml_h::ensure_valid_string_option(opt_val, options);
}

template <typename Key = std::string>
bool contains(toml::value const &v, Key const &k) {
  return v.contains(k);
}

template <typename Key = std::string>
bool contains_rec(toml::value const &v, std::vector<Key> const &ks) {
  if (!v.contains(ks.front())) {
    return false;
  }

  if (ks.size() > 1) {
    return toml_h::contains_rec(v.at(ks.front()),
                                std::vector<Key>(ks.begin() + 1, ks.end()));
  } else {
    return v.contains(ks.front());
  }
}

} // namespace toml_h
