#include "gtest/gtest.h"

#include "grid.hpp"

TEST(Grid, CheckPeriodicMappingOfCoordinates) {
    EXPECT_EQ(artemis::grid::periodic_mapping(-6, 6), 0);
    EXPECT_EQ(artemis::grid::periodic_mapping(-5, 6), 1);
    EXPECT_EQ(artemis::grid::periodic_mapping(-4, 6), 2);
    EXPECT_EQ(artemis::grid::periodic_mapping(-3, 6), 3);
    EXPECT_EQ(artemis::grid::periodic_mapping(-2, 6), 4);
    EXPECT_EQ(artemis::grid::periodic_mapping(-1, 6), 5);
    EXPECT_EQ(artemis::grid::periodic_mapping(0, 6), 0);
    EXPECT_EQ(artemis::grid::periodic_mapping(1, 6), 1);
    EXPECT_EQ(artemis::grid::periodic_mapping(2, 6), 2);
    EXPECT_EQ(artemis::grid::periodic_mapping(3, 6), 3);
    EXPECT_EQ(artemis::grid::periodic_mapping(4, 6), 4);
    EXPECT_EQ(artemis::grid::periodic_mapping(5, 6), 5);
    EXPECT_EQ(artemis::grid::periodic_mapping(6, 6), 0);
    EXPECT_EQ(artemis::grid::periodic_mapping(7, 6), 1);
    EXPECT_EQ(artemis::grid::periodic_mapping(8, 6), 2);
}

TEST(Grid, CellConversionOfCoordinatesToCellIds) {
    EXPECT_EQ(artemis::grid::convert_to_cell_id({-3, 0, 0}, {4, 3, 2}), 1);
    EXPECT_EQ(artemis::grid::convert_to_cell_id({-2, 0, 0}, {4, 3, 2}), 2);
    EXPECT_EQ(artemis::grid::convert_to_cell_id({-1, 0, 0}, {4, 3, 2}), 3);
    EXPECT_EQ(artemis::grid::convert_to_cell_id({0, 0, 0}, {4, 3, 2}), 0);
    EXPECT_EQ(artemis::grid::convert_to_cell_id({1, 0, 0}, {4, 3, 2}), 1);
    EXPECT_EQ(artemis::grid::convert_to_cell_id({2, 0, 0}, {4, 3, 2}), 2);
    EXPECT_EQ(artemis::grid::convert_to_cell_id({3, 0, 0}, {4, 3, 2}), 3);
    EXPECT_EQ(artemis::grid::convert_to_cell_id({0, 1, 0}, {4, 3, 2}), 4);
    EXPECT_EQ(artemis::grid::convert_to_cell_id({0, 0, 1}, {4, 3, 2}), 12);
}

TEST(Grid, ConversionOfCellIdsToCoordinates) {

    // EXPECT_EQ(artemis::grid::convert_to_coordinates(0, {4, 3, 2}), {0, 0, 0});
}