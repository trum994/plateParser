import os
import sys
import pandas as pd


class WellPlate96:

    def __init__(self, well_count=96):
        self.well_count = well_count
        # creating 96 column names in a 8x12 nested loop
        columns = "ABCDEFGH"
        rfu_columns = []
        for char in columns:
            for i in range(1, 13):
                column_name = char + str(i)
                rfu_columns.append(column_name)
        self.this_df = pd.DataFrame(columns=rfu_columns)
        self.labels = self.this_df.iloc[:, 0]

    def get_exp_count(self):
        return self.this_df.shape[0]

    def load_raw_data(self, raw_in, this_is_folder):
        # read in the file skipping header including column names and footer rows
        some_df = pd.read_csv(raw_in, skiprows=3, skipfooter=2, header=None, sep='\t', engine='python',
                              encoding_errors='ignore')

        # file format creates two empty columns at the end, drop these
        some_df = some_df.iloc[:, :-2]

        # If input is coming from folder time is not in file, so parse from filename and inject
        if this_is_folder:
            time_from_filename = os.path.splitext(os.path.basename(raw_in))[0].replace("_", ".")
            self.labels = self.labels.append(pd.Series([time_from_filename]), ignore_index=True)
        else:
            # get the time intervals
            just_times = some_df.iloc[:, 0]
            just_times.dropna(axis=0, how="any", inplace=True)
            self.labels = self.labels.append(pd.Series(just_times))

        # drop first two columns and we now only have pure data, ie. relative fluorescence units RFU
        some_df = some_df.iloc[:, 2:]

        # scroll in increments of 9 to account for blank row at the end of each chunk
        for i in range(0, some_df.shape[0], 9):
            # only take top 8 rows, 9th row is blank
            chunk = some_df.iloc[i:i+8, :]
            my_series = pd.Series(chunk.to_numpy(copy=True).flatten(), index=self.this_df.columns)
            self.this_df = self.this_df.append(my_series, ignore_index=True)

    def get_columns(self, my_coord_list):
        print(self.this_df[my_coord_list].to_csv(float_format='%.3f'))


class CoordList:

    def __init__(self, user_in):
        self.user_in = user_in
        self.fixed_list = []
        self.coordinates_to_cols()

    def coordinates_to_cols(self):
        final_list = []
        for item in self.user_in.split(","):
            if "-" in item:
                # process further since it contains a dash
                final_list = final_list + self.dash_item_to_cols(item)
            else:
                # add to the list as is
                final_list.append(item)
        self.fixed_list = final_list

    @staticmethod
    def dash_item_to_cols(dashed_item):
        dashed_list = []
        first_part = dashed_item.split("-")[0]
        second_part = dashed_item.split("-")[1]

        if (first_part.isalpha()) and (len(first_part) == 1):
            # this means we loop around characters
            begin = first_part
            end = second_part[0]
            constant = second_part[1:]
            for my_char_num in range(ord(begin), ord(end)+1):
                column_name = chr(my_char_num) + constant
                dashed_list.append(column_name)
        else:
            # this means we loop around numbers
            constant = first_part[0]
            begin = int(first_part[1:])
            end = int(second_part)
            for my_number in range(begin, end+1):
                column_name = constant + str(my_number)
                dashed_list.append(column_name)
        return dashed_list


def run_parser(my_df, my_coord):
    print("Running parser")
    print("Will process with coordinates: " + my_coord)
    # Create coordinates object
    my_co = CoordList(my_coord)
    print("Fixed list of coordinates: {0}".format(str(my_co.fixed_list)))

    # Print number of time points
    print("Read total of " + str(my_df.get_exp_count()) + " time points.")

    # Add times as an index for the data frame
    my_df.labels.name = "Time"
    my_df.this_df.index = my_df.labels
    my_df.get_columns(my_co.fixed_list)


def run_cytotox(my_df):
    print("Running cytotox")

    bgr_wells = "B-D1"
    d1c_wells = "B-D2"
    d2c_wells = "E-G2"
    d1e_wells = "B2-11,C2-11,D2-11"
    d2e_wells = "E2-11,F2-11,G2-11"

    # Create coordinates for background, drug1, drug2 wells
    bgr_co = CoordList(bgr_wells)
    d1c_co = CoordList(d1c_wells)
    d2c_co = CoordList(d2c_wells)
    d1e_co = CoordList(d1e_wells)
    d2e_co = CoordList(d2e_wells)

    # print("Fixed list of bgr coordinates: {0}".format(str(bgr_co.fixed_list)))
    # print("Fixed list of d1c coordinates: {0}".format(str(d1c_co.fixed_list)))
    # print("Fixed list of d2c coordinates: {0}".format(str(d2c_co.fixed_list)))
    # print("Fixed list of d1e coordinates: {0}".format(str(d1e_co.fixed_list)))
    # print("Fixed list of d2e coordinates: {0}".format(str(d2e_co.fixed_list)))

    my_df.labels.name = "Cell-Drug"
    my_df.this_df.index = my_df.labels

    # Calculate average background and add as new column at the end
    my_df.this_df['BG'] = my_df.this_df[bgr_co.fixed_list].mean(axis=1)

    # Subtract background from all values
    my_df.this_df = my_df.this_df.sub(my_df.this_df['BG'], axis=0)

    # Calculate drug1 and drug2 control values and add as new column at the end
    my_df.this_df['DR1C'] = my_df.this_df[d1c_co.fixed_list].mean(axis=1)
    my_df.this_df['DR2C'] = my_df.this_df[d2c_co.fixed_list].mean(axis=1)

    # Calculate percentage by dividing everything by control values and multiplying by 100
    my_df.this_df[d1e_co.fixed_list] = my_df.this_df[d1e_co.fixed_list].div(my_df.this_df['DR1C'], axis=0).multiply(100)
    my_df.this_df[d2e_co.fixed_list] = my_df.this_df[d2e_co.fixed_list].div(my_df.this_df['DR2C'], axis=0).multiply(100)

    final_df = pd.DataFrame()
    # Each row is an input file. Now converting each file back into a table/dataframe
    for index, row in my_df.this_df.iterrows():
        drug1_array = row[d1e_co.fixed_list].to_frame().transpose().values.reshape(3, 10)
        drug1_df = pd.DataFrame(drug1_array, columns=range(2, 12))
        drug1_df.index = list('BCD')
        drug1_df['drug'] = "drug1"
        drug1_df['filename'] = index
        final_df = final_df.append(drug1_df)

        drug2_array = row[d2e_co.fixed_list].to_frame().transpose().values.reshape(3, 10)
        drug2_df = pd.DataFrame(drug2_array, columns=range(2, 12))
        drug2_df.index = list('EFG')
        drug2_df['drug'] = "drug2"
        drug2_df['filename'] = index
        final_df = final_df.append(drug2_df)

    print(final_df)
    final_df.to_excel("output.xlsx")


def get_file_list(my_input_list):
    this_file_list = []
    if os.path.isdir(my_input_list):
        for filename in os.listdir(my_input_list):
            if filename.endswith(".txt"):
                filepath = os.path.join(my_input_list, filename)
                this_file_list.append(filepath)
    else:
        this_file_list.append(input_list)
    print(sorted(this_file_list))
    return sorted(this_file_list)


if __name__ == '__main__':

    if len(sys.argv) != 3:
        print("Required both arguments: input file(s) and comma separated coordinates or just cytotox")
        sys.exit(1)
    input_list = sys.argv[1]

    file_list = get_file_list(input_list)
    print("Will process file(s): " + input_list)

    # Initialize 96 well plate with columns
    the_df = WellPlate96()

    for this_file in sorted(file_list):
        # Load raw data in dataframe
        the_df.load_raw_data(this_file, os.path.isdir(input_list))

    mode = ""
    if sys.argv[2] == "cytotox":
        mode = "cytotox"
        run_cytotox(the_df)
    else:
        mode = "parser"
        run_parser(the_df, my_coord=sys.argv[2])
