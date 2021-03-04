#ifndef STRING_UTILS_H
#define STRING_UTILS_H 1

#include <string>
#include <vector>
#include <sstream>
#include <algorithm>

static inline std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

static inline std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
}

// trim from start
static inline std::string &ltrim(std::string &my_string) {
        my_string.erase(my_string.begin(), std::find_if(my_string.begin(), my_string.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
        return my_string;
}

// trim from end
static inline std::string &rtrim(std::string &my_string) {
        my_string.erase(std::find_if(my_string.rbegin(), my_string.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), my_string.end());
        return my_string;
}

// trim from both ends
static inline std::string &trim(std::string &my_string) {
        return ltrim(rtrim(my_string));
}


static inline std::string remove_extension(const std::string &filename) {
    std::vector<std::string> temp_vector;
    std::string final_string = "";
    temp_vector = split(filename, '.');
    final_string = temp_vector[0];
    // The extension is the last token of a filename tokenized by "."
    // so add all tokens together except for the last one
    for (size_t i = 1; i < temp_vector.size() - 1; i++) {
        final_string += "." + temp_vector[i];
    }

    return final_string;
}

static inline std::string find_extension(const std::string &filename) {
    std::vector<std::string> temp_vector;
    std::string final_string = "";
    temp_vector = split(filename, '.');
    if (temp_vector.size() > 0) {
        return temp_vector[temp_vector.size() - 1];
    } else {
        return "";
    }
}


#endif
