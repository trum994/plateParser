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


if __name__ == '__main__':

    if len(sys.argv) != 3:
        print("Required both arguments: input file(s) and comma separated coordinates")
        sys.exit(1)
    input_list = sys.argv[1]
    coordinates = sys.argv[2]
    bg_wells = "B-D1"
    d1_wells = "B-D2"
    d2_wells = "E-G2"

    file_list = []
    print("Will process file(s): " + input_list + " with coordinates: " + coordinates)

    is_folder = False
    if os.path.isdir(input_list):
        is_folder = True
        for filename in os.listdir(input_list):
            if filename.endswith(".txt"):
                filepath = os.path.join(input_list, filename)
                file_list.append(filepath)
    else:
        file_list.append(input_list)
    print(sorted(file_list))

    # Initialize 96 well plate with columns
    my_df = WellPlate96()

    # Create coordinates object
    my_co = CoordList(coordinates)
    print("Fixed list of coordinates: {0}".format(str(my_co.fixed_list)))

    # Create coordinates for background, drug1, drug2 wells
    bg_co = CoordList(bg_wells)
    d1_co = CoordList(d1_wells)
    d2_co = CoordList(d2_wells)
    print("Fixed list of bg coordinates: {0}".format(str(bg_co.fixed_list)))
    print("Fixed list of d1 coordinates: {0}".format(str(d1_co.fixed_list)))
    print("Fixed list of d2 coordinates: {0}".format(str(d2_co.fixed_list)))

    for this_file in sorted(file_list):
        # Load raw data in dataframe
        my_df.load_raw_data(this_file, is_folder)

    # Print number of time points
    print("Read total of " + str(my_df.get_exp_count()) + " time points.")

    # Add times as an index for the data frame
    my_df.labels.name = "Time"
    my_df.this_df.index = my_df.labels
    my_df.get_columns(my_co.fixed_list)

    # Calculate average background and add as new column at the end
    my_df.this_df['BG'] = my_df.this_df[bg_co.fixed_list].mean(axis=1)
    print(my_df.this_df)

    # Subtract background from all values
    my_df.this_df = my_df.this_df.sub(my_df.this_df['BG'], axis=0)
    print(my_df.this_df)

    # Calculate drug1 and drug2 averages and as new column at the end
    my_df.this_df['DR1'] = my_df.this_df[d1_co.fixed_list].mean(axis=1)
    my_df.this_df['DR2'] = my_df.this_df[d2_co.fixed_list].mean(axis=1)
    print(my_df.this_df)
    print(my_df.this_df[['B2', 'C2', 'D2']])
