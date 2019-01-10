#pragma once

#include <string>

class ReadItem {
    public:
        std::string sequence1;
        std::string sequence2;
        bool paired = false;
        ReadItem(ReadItem const&) = delete;
        void operator=(ReadItem const&) = delete;
        ReadItem(const std::string &);
        ReadItem(const std::string &, const std::string &);
};
