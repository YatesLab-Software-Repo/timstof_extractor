import timsdata, sqlite3, sys, time, os
import numpy as np, matplotlib.pyplot as plt
from scipy import constants
import sqlite3
from sqlite3 import Error
from pathlib import Path
import glob
import time

def K0toCCS (K0, q, m_ion, m_gas, T):
    mu = m_ion*m_gas/(m_ion+m_gas)
    T0 = 273.15
    p0 = 1.0132e5 # 1 atm
    N0 = p0 / (constants.k * T0)
    return (3.0/16.0) * (1/N0) * np.sqrt(2 * np.pi/ (mu*constants.k * T)) * (q/K0)


def enhance_signal(intensity):
    return (intensity**1.414+100.232)**1.414

def create_connection(analysis_dir):
    import os
    db_file = os.path.join(analysis_dir, 'analysis.tdf')
    """ create a database connection to the SQLite database
        specified by the db_file
    :param db_file: database file
    :return: Connection object or None
    """
    try:
        conn = sqlite3.connect(db_file)
        return conn
    except Error as e:
        print(e)

    return None


def select_all_PasefFrameMsMsInfo(conn):
    """
    Query all rows in the tasks table
    :param conn: the Connection object
    :return:
    """
    cur = conn.cursor()
    cur.execute("SELECT * FROM PasefFrameMsMsInfo")

    rows = cur.fetchall()

    data_all = np.array(rows)
    return data_all


def select_all_Frames(conn):
    """
    Query all rows in the tasks table
    :param conn: the Connection object
    :return:
    """
    cur = conn.cursor()
    cur.execute("SELECT * FROM Frames")

    rows = cur.fetchall()

    data_all = np.array(rows)
    return data_all


def select_all_Precursors(conn):
    """
    Query all rows in the tasks table
    :param conn: the Connection object
    :return:
    """
    cur = conn.cursor()
    cur.execute("SELECT * FROM Precursors")

    rows = cur.fetchall()

    data_all = np.array(rows)
    return data_all


def msms_frame_parent_dict(all_frame):
    i = 1
    parent_frame_id_dict = {}
    msms_type_array = np.array(all_frame).transpose()[4]
    last_zero = 0
    for each in msms_type_array:
        if each == '0':
            parent_frame_id_dict[i] = i
            last_zero = i
        elif each == '8':
            parent_frame_id_dict[i] = last_zero
        i += 1
    return parent_frame_id_dict

def runTimstofConversiont(input, output=''):
    analysis_dir = input

    td = timsdata.TimsData(analysis_dir)
    conn = td.conn

    # create a database connection
    conn = create_connection(analysis_dir)

    with conn:
        # print("2. Query all tasks")
        msms_data = select_all_PasefFrameMsMsInfo(conn)
        # print msms_data[0:5]
        all_frame = select_all_Frames(conn)
        # print all_frame[0:5]
        precursor_list = select_all_Precursors(conn)

    all_ms1_frames = [a for a in all_frame if a[4] == '0']

    precursor_array = np.array(precursor_list)  # 'ID', 'LargestPeakMz', 'AverageMz', 'MonoisotopicMz', 'Charge', 'ScanNumber', 'Intensity', 'Parent'
    frame_parent_dict = msms_frame_parent_dict(all_frame)

    # precursor_df = pd.DataFrame(data=precursor_array, index=precursor_array[:, 0], columns=['ID', 'LargestPeakMz', 'AverageMz', 'MonoisotopicMz', 'Charge', 'ScanNumber', 'Intensity', 'Parent'])

    parent_frame_array = np.array(precursor_array[:, 7])
    skip_ms2 = False
    frame_index_list = []
    last_val = 0
    for idx, val in enumerate(parent_frame_array):
        if val != last_val:
            frame_index_list.append(idx)
        last_val = val
    frame_index_list.append(idx + 1)
    frame_start_end_dict = {}

    for idx, val in enumerate(frame_index_list[:-1]):
        frame_start_end_dict[parent_frame_array[val]] = (frame_index_list[idx], frame_index_list[idx + 1])

    ms2_header = 'H\tExtractor\tTimsTOF_extractor\nH\tExtractorVersion\t0.0.1\nH\tComments\tTimsTOF_extractor written by Yu Gao, 2018\nH\tExtractorOptions\tMSn\nH\tAcquisitionMethod\tData-Dependent\nH\tInstrumentType\tTIMSTOF\nH\tDataType\tCentroid\nH\tScanType\tMS2\nH\tResolution\nH\tIsolationWindow\nH\tFirstScan\t1\nH\tLastScan\t%s\nH\tMonoIsotopic PrecMz\tTrue\n' % len(
        msms_data)
    ms2_file_name=os.path.basename(analysis_dir).split('.')[0]+'_nopd.ms2'
    ms1_file_name = os.path.basename(analysis_dir).split('.')[0]+'_nopd.ms1'
    ms1_scan_set = set()
    if len(output)> 0:
        ms2_file_name = output
        ms1_file_name = output.replace('.ms2', '.ms1')
    else:
        os.chdir(sys.argv[2])
    if not skip_ms2:
        with open(ms2_file_name, 'w') as output_file:
            output_file.write(ms2_header)
            progress = 0
            last_parent_frame = 0
            for each_msms_data in msms_data:
                frame_id, scan_begin, scan_end, precursor_mass, isolation_width, collision_energy, precursor_seq = each_msms_data
                ms_ms = td.readPasefMsMs([precursor_seq]);

                # print frame_id,scan_begin, scan_end, precursor_mass
                frame_id_int = int(frame_id)
                scan_begin_int = int(scan_begin)
                scan_end_int = int(scan_end)
                precursor_seq_int = int(precursor_seq)
                scans = td.readScans(frame_id_int, scan_begin_int, scan_end_int)
                index_intensity = np.concatenate(scans, axis=1)






               # print("scans: ", index_intensity);
                mass_array = td.indexToMz(frame_id_int, index_intensity[0])
                temp = np.array(list(zip(mass_array, index_intensity[1])))
                mass_intensity = np.around(temp, decimals=4)
                if len(mass_intensity) == 0:
                    continue
                sorted_mass_intensity = mass_intensity[mass_intensity[:, 0].argsort()]
                # mass_intensity.sort(key=lambda x: x[0])

                mass_charge_list = precursor_array[int(each_msms_data[6] - 1), [1, 2, 3, 4, 5]]

                if mass_charge_list[2] == None:
                    mass_charge_list[2] = mass_charge_list[0]
                if mass_charge_list[3] == None:
                    mass_charge_list[3] = 0

                scan_no = int(mass_charge_list[4])

                if mass_charge_list[3] != 0:
                    mz = (mass_charge_list[2] - 1.0072766) * mass_charge_list[3] + 1.0072766
                    output_file.write("S\t%06d\t%06d\t%.4f\n" % (progress, progress, mass_charge_list[2]))
                    output_file.write("I\tRetTime\t%.2f\n" % float(all_frame[frame_id_int - 1][1]))
                    output_file.write("Z\t%d\t%.4f\n" % (mass_charge_list[3], mz))





                else:
                    mz_charge2 = (mass_charge_list[2] - 1.0072766) * 2 + 1.0072766
                    mz_charge3 = (mass_charge_list[2] - 1.0072766) * 3 + 1.0072766
                    output_file.write("S\t%06d\t%06d\t%.4f\n" % (progress, progress, mass_charge_list[2]))
                    output_file.write("I\tRetTime\t%.2f\n" % float(all_frame[frame_id_int - 1][1]))
                    output_file.write("Z\t%d\t%.4f\n" % (2, mz_charge2))
                    output_file.write("Z\t%d\t%.4f\n" % (3, mz_charge3))

                last_mass = sorted_mass_intensity[0][0]
                total_intensity = 0

                for i, each_pair in enumerate(sorted_mass_intensity):
                    mass, intensity = each_pair
                    if mass == last_mass:
                        total_intensity += intensity
                    else:
                        output_file.write("%.4f %.1f \n" % (last_mass, enhance_signal(total_intensity)))
                        total_intensity = intensity
                    last_mass = mass

                output_file.write("%.4f %.1f\n" % (last_mass, enhance_signal(total_intensity)))


                progress += 1
                if progress % 5000 == 0:
                    print("progress ms2: %.1f%%" % (float(progress) / len(msms_data) * 100), time.clock() - start_time)

    with open(ms1_file_name, 'w') as output_file:
        output_file.write(ms2_header)
        progress = 0
        for i, frame in enumerate(all_ms1_frames):
            id = int(frame[0])
            num_scans = int(frame[8])
            #time = frame[1]
            index_intensity_arr = td.readScans(id, 0, num_scans)
            index_intensity_carr = np.concatenate(index_intensity_arr, axis=1)
            mobility_index = [i for i, row in enumerate(index_intensity_arr) for j in range(len(row[0]))]

            mass_array = td.indexToMz(id, index_intensity_carr[0])
            one_over_k0 = td.scanNumToOneOverK0(id, mobility_index)

            temp = np.array(list(zip(mass_array, index_intensity_carr[1], one_over_k0)))
            mass_intensity = np.around(temp, decimals=4)
            sorted_mass_intensity = mass_intensity[mass_intensity[:, 0].argsort()]

            rt_time = 0 if i == 0 else all_frame[i-1][1]

            output_file.write("S\t%06d\t%06d\n" % (id, id))
            output_file.write("I\tRetTime\t%.2f\n" % float(rt_time))

            for i, row in enumerate(sorted_mass_intensity):
                output_file.write("%.4f %.1f %.4f\n" % (row[0], row[1],
                                                        row[-1]))
            progress += 1
            if progress % 5000 == 0:
                print("progress ms1 %.1f%%" % (float(progress) / len(all_frame) * 100), time.clock() - start_time)



if __name__ == '__main__':
    if len(sys.argv)==1:
        print("Usage: extract_msn_nopd [source data directory (.d)] [target directory for output]")
    else:
        analysis_dir = sys.argv[1]

    start_time = time.clock()

    dirs_to_analyze = []
    if analysis_dir[-1] == '*':
        #p = Path('.')
        #print(p)
        for f in glob.glob(analysis_dir):
            if os.path.isdir(f):
                tdf_file = os.path.join(f, "analysis.tdf")
                if os.path.exists(tdf_file):
                    print(tdf_file)
                    dirs_to_analyze.append(f)

        for input in dirs_to_analyze:
            ms2_file_name = os.path.basename(input).split('.')[0] + '_nopd.ms2'
            ms2_file_name = os.path.join(input,ms2_file_name)
            print(ms2_file_name)
            output = input + os.path.sep + ms2_file_name
            runTimstofConversiont(input,ms2_file_name)
    else:
        runTimstofConversiont(analysis_dir)


