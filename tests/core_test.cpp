#include <gtest/gtest.h>
#include "mphhc/core.hpp"

TEST(CoreTest, AddFunction) {
    EXPECT_EQ(mphhc::add(1, 2), 3);
    EXPECT_NE(mphhc::add(1, 2), 4);
}

TEST(CoreTest, VersionString) {
    EXPECT_EQ(mphhc::get_version(), "0.1.0");
}
