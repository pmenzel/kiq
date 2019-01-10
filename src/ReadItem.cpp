#include "ReadItem.hpp"

ReadItem::ReadItem(const std::string & s1) : sequence1(s1) { }

ReadItem::ReadItem(const std::string & s1, const std::string & s2) : sequence1(s1), sequence2(s2), paired(true) { }

