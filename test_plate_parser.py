import pytest
from plate_parser import WellPlate96, CoordList


@pytest.fixture
def my_fix_plate():
    my_test_plate = WellPlate96()
    return my_test_plate


def test_well_count(my_fix_plate):
    assert my_fix_plate.well_count == 96


def test_col_names(my_fix_plate):
    assert (str(my_fix_plate.this_df.columns.shape) == "(96,)")


def test_load_raw_data(my_fix_plate):
    raw_in = "rawInput.txt"
    my_fix_plate.load_raw_data(raw_in, False)
    assert my_fix_plate.this_df.iloc[0, 0] == 40.144
    assert my_fix_plate.this_df.iloc[-1, -1] == 37.01


@pytest.fixture()
def my_fix_coordinates_with_dash():
    my_test_coordinates_with_dash = CoordList("A2-3,A-D4")
    return my_test_coordinates_with_dash


@pytest.fixture()
def my_fix_coordinates_standard():
    my_test_coordinates_standard = CoordList("F3,G7,B5")
    return my_test_coordinates_standard


def test_coordinates_to_cols_with_dash(my_fix_coordinates_with_dash):
    assert my_fix_coordinates_with_dash.fixed_list[0] == "A2"
    assert my_fix_coordinates_with_dash.fixed_list[-1] == "D4"


def test_coordinates_to_cols_standard(my_fix_coordinates_standard):
    assert my_fix_coordinates_standard.fixed_list[0] == "F3"
    assert my_fix_coordinates_standard.fixed_list[-1] == "B5"
