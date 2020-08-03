import timsdata, sqlite3, sys, time, os
import numpy as np, matplotlib.pyplot as plt
from scipy import constants
import sqlite3
from sqlite3 import Error
import glob
import time
import math

place_high = 3
precursor_counter = 0
convert_ms2 = True
convert_ms1 = True
version = "0.1.2"


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


def generate_MS1_scan(conn, id, get_last =False):
    global place_high
    cur = conn.cursor();
    key = (id,)
    cur.execute('''
        select Frame, ScanNumBegin from PasefFrameMsMsInfo 
        INNER JOIN Precursors on PasefFrameMsMsInfo.Precursor = Precursors.Id 
        where Precursors.Parent = ? order by frame, ScanNumBegin ''', key)
    rows = cur.fetchall()
    if len(rows) > 0:
        index = len(rows) - 1 if get_last else 0
        data = np.array(rows[index])
        tims_scan_num_begin = int(data[1])
        frame = int(data[0])
        place = math.ceil(math.log10(tims_scan_num_begin))
        if place > place_high:
            place_high = place
        else:
            place = place_high

        ms1_scan = frame * 10 ** place + tims_scan_num_begin-1
        return ms1_scan
    else:
        return -1


def generate_ms1_scan_v2(frame_id, precursor_list):
    global precursor_counter
    while precursor_counter < len(precursor_list) and precursor_list[precursor_counter][-1] < frame_id:
        precursor_counter += 1
    if precursor_counter >= len(precursor_list):
        id = precursor_list[len(precursor_list)-1][0]
        parent = precursor_list[len(precursor_list)-1][-1]
        return id+parent+1, id, parent
    elif precursor_list[precursor_counter][-1] > frame_id and precursor_counter > 0:
        id = precursor_list[precursor_counter-1][0]
        parent = precursor_list[precursor_counter-1][-1]
        return id+parent+1, id, parent
    else:
        id = precursor_list[precursor_counter][0]
        parent = precursor_list[precursor_counter][-1]
        return id+parent-1,  id, parent


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


def build_offset_map(precursor_map, all_ms1_list):
    offset = 0
    offset_map = {}
    for row in all_ms1_list:
        frame_id = int(row[0])
        offset_map[frame_id] = offset
        if frame_id in precursor_map:
            precursor_id = int(precursor_map[frame_id][0])
            parent = int(precursor_map[frame_id][-1])
            offset += precursor_id + parent
    return offset_map


def build_frame_id_ms1_scan_map(precursor_map, all_ms1_list):
    frame_id_ms1_scan_map = {}
    ms2_map = {}
    prev_scan = 0
    for row in all_ms1_list:
        frame_id = int(row[0])
        prev_scan += 1
        frame_id_ms1_scan_map[frame_id] = prev_scan
        if frame_id in precursor_map:
            if frame_id not in ms2_map:
                ms2_map[frame_id] = {}
            for count, rows in enumerate(precursor_map[frame_id], prev_scan+1):
                prec_id = int(rows[0])
                ms2_map[frame_id][prec_id] = count
            prev_scan += len(precursor_map[frame_id])
    return frame_id_ms1_scan_map, ms2_map


def run_timstof_conversion(input, output=''):
    global place_high
    global precursor_counter
    analysis_dir = input

    td = timsdata.TimsData(analysis_dir)
    conn = td.conn

    # create a database connection
    #conn = create_connection(analysis_dir)
    precursor_map ={}

    with conn:
        # print("2. Query all tasks")
        msms_data = select_all_PasefFrameMsMsInfo(conn)
        # print msms_data[0:5]
        all_frame = select_all_Frames(conn)
        # print all_frame[0:5]
        precursor_list = select_all_Precursors(conn)
        for row in precursor_list:
            parent_id = int(row[-1])
            if parent_id not in precursor_map:
                precursor_map[parent_id] = []
            precursor_map[parent_id].append(row)


    all_ms1_frames = [a for a in all_frame if a[4] == '0']

    frame_id_ms1_scan_map, ms2_scan_map = build_frame_id_ms1_scan_map(precursor_map, all_ms1_frames)
    #offset_map = build_offset_map(precursor_map, all_ms1_frames)

    precursor_array = np.array(precursor_list)   # 'ID', 'LargestPeakMz', 'AverageMz', 'MonoisotopicMz', 'Charge', 'ScanNumber', 'Intensity', 'Parent'
    # frame_parent_dict = msms_frame_parent_dict(all_frame)
    # parent_ms2_scan_map = build_frame_id_last_ms2_scan_map(precursor_list)



    parent_frame_array = np.array(precursor_array[:, 7])
    frame_index_list = []
    last_val = 0
    for idx, val in enumerate(parent_frame_array):
        if val != last_val:
            frame_index_list.append(idx)
        last_val = val
#    frame_index_list.append(idx + 1)
  #  frame_start_end_dict = {}

    #for idx, val in enumerate(frame_index_list[:-1]):
       # frame_start_end_dict[parent_frame_array[val]] = (frame_index_list[idx], frame_index_list[idx + 1])

    ms2_header = 'H\tExtractor\tTimsTOF_extractor\n' \
                 'H\tExtractorVersion\t{}\n' \
                 'H\tPublicationDate\t20-02-2020\n' \
                 'H\tComments\tTimsTOF_extractor written by Yu Gao, 2018\n' \
                 'H\tComments\tTimsTOF_extractor modified by Titus Jung, 2019\n' \
                 'H\tExtractorOptions\tMSn\n' \
                 'H\tAcquisitionMethod\tData-Dependent\n' \
                 'H\tInstrumentType\tTIMSTOF\n' \
                 'H\tDataType\tCentroid\n' \
                 'H\tScanType\tMS2\n' \
                 'H\tResolution\n' \
                 'H\tIsolationWindow\n' \
                 'H\tFirstScan\t1\n' \
                 'H\tLastScan\t{}\n' \
                 'H\tMonoIsotopic PrecMz\tTrue\n'.format(version, len(msms_data))

    ms1_header = 'H\tExtractor\tTimsTOF_extractor\n' \
                 'H\tExtractorVersion\t{}\n' \
                 'H\tPublicationDate\t20-02-2020\n' \
                 'H\tComments\tTimsTOF_extractor written by Yu Gao, 2018\n' \
                 'H\tComments\tTimsTOF_extractor modified by Titus Jung, 2019\n' \
                 'H\tExtractorOptions\tMSn\n' \
                 'H\tAcquisitionMethod\tData-Dependent\n' \
                 'H\tInstrumentType\tTIMSTOF\n' \
                 'H\tScanType\tMS1\n' .format(version)

    ms2_file_name = os.path.basename(analysis_dir).split('.')[0]+'_nopd.ms2'
    ms1_file_name = os.path.basename(analysis_dir).split('.')[0]+'_nopd.ms1'
    ms1_scan_set = set()
    if len(output) > 0:
        ms2_file_name = output
        ms1_file_name = output.replace('.ms2', '.ms1')
    #else:
       # os.chdir(sys.argv[2])
    if convert_ms2:
        with open(ms2_file_name, 'w') as output_file:
            output_file.write(ms2_header)
            progress = 0
            for row in precursor_list:
                prc_id, largest_preak_mz, average_mz, monoisotopic_mz, cs, scan_number, intensity, parent = row
                prc_id_int = int(prc_id)
                if monoisotopic_mz is not None and cs is not None:
                    prc_mass_mz = float(monoisotopic_mz)
                    prc_mass = (prc_mass_mz*cs)-(cs-1)*1.007276466

                    mz_int_arr = td.readPasefMsMs([prc_id_int])
                    parent_index = int(parent)
                    scan_id = ms2_scan_map[parent_index][prc_id_int]
                    rt_time = float(all_frame[parent_index][1])
                    k0 = td.scanNumToOneOverK0(parent_index, [scan_number])
                    mz_arr = mz_int_arr[prc_id_int][0]
                    if len(mz_arr) > 0:
                        output_file.write("S\t{0:06d}\t{1:06d}\t{2:.4f}\n".format(scan_id, scan_id, prc_mass_mz))
                        output_file.write("I\tTIMSTOF_Parent_ID\t{}\n".format(parent))
                        output_file.write("I\tTIMSTOF_Precursor_ID\t{}\n".format(prc_id))
                        output_file.write("I\tRetTime\t{0:.4f}\n".format(rt_time))
                        output_file.write("I\tIon Mobility\t{0:.4f}\n".format(k0[0]))
                        output_file.write("Z\t{1}\t{0:.4f}\n".format(prc_mass, cs))

                        int_arr = mz_int_arr[prc_id_int][1]
                        for j in range(0, len(mz_arr)):
                            output_file.write("%.4f %.1f \n" % (mz_arr[j], int_arr[j]))

                    progress += 1
                    if progress % 5000 == 0:
                        print("progress ms2: %.1f%%" % (float(progress) / len(precursor_list) * 100), time.process_time() - start_time)
    if convert_ms1:
        with open(ms1_file_name, 'w') as output_file:
            output_file.write(ms1_header)
            progress = 0
            prev_id = 0
            #scan_set = set()
            prev_scan = 0
            precursor_counter = 0
            lines = []
            for i, frame in enumerate(all_ms1_frames):
                id = int(frame[0])
                num_scans = int(frame[8])

                index_intensity_arr = td.readScans(id, 0, num_scans)
                index_intensity_carr = np.concatenate(index_intensity_arr, axis=1)
                mobility_index = [i for i, row in enumerate(index_intensity_arr) for j in range(len(row[0]))]

                mass_array = td.indexToMz(id, index_intensity_carr[0])
                one_over_k0 = td.scanNumToOneOverK0(id, mobility_index)
                voltage = td.scanNumToVoltage(id, mobility_index)
                temp = np.array(list(zip(mass_array, index_intensity_carr[1], one_over_k0, voltage)))
                mass_intensity = np.around(temp, decimals=4)
                sorted_mass_intensity = mass_intensity[mass_intensity[:, 0].argsort()]
                scan_num = frame_id_ms1_scan_map[id]
                if len(sorted_mass_intensity) >0:
                    rt_time = 0 if i == 0 else all_ms1_frames[i-1][1]
                    lines.append("S\t%06d\t%06d\n" % (scan_num, scan_num))
                    lines.append("I\tTIMSTOF_Frame_id\t{}\n".format(id))
                    lines.append("I\tRetTime\t%.2f\n" % float(rt_time))
                    for row in sorted_mass_intensity:
                        x_str = "%.4f %.1f %.4f \n" % (row[0], row[1],  row[-2])
                        lines.append(x_str)
                # output_file.write("S\t%06d\t%06d\n" % (scan_num, scan_num))
                # output_file.write("I\tTIMSTOF_Frame_id\t{}\n".format(id))
                # output_file.write("I\tRetTime\t%.2f\n" % float(rt_time))
                # output_file.writelines("%.4f %.1f %.4f\n" % (row[0], row[1],
                # row[-1]) for row in sorted_mass_intensity)
                if len(lines) > 1_000_000:
                    output_file.writelines(lines)
                    lines = []

                progress += 1
                if progress % 5000 == 0:
                    print("progress ms1 %.1f%%" % (float(progress) / len(all_ms1_frames) * 100), time.process_time()
                          - start_time)
            output_file.writelines(lines)
            lines = []
    conn.close()

def is_valid_timstof_dir(path):
    if os.path.isdir(path):
        tdf_file = os.path.join(path, "analysis.tdf")
        if os.path.exists(tdf_file):
            return True
    return False


if __name__ == '__main__':
    dirs_to_analyze = []
    print("Running timSTOFExtractor {}".format(version))
    if len(sys.argv) == 1:
        print("Usage: extract_msn_nopd [source data directory (.d)]")
        print("--skip-ms1: skips creating ms1 files")
        print("--skip-ms2: skips creating ms2 files")
    else:
        analysis_dir = sys.argv[1]
    for i, x in enumerate(sys.argv):
        if x == "--skip-ms1":
            print("<<<<: skipping ms1")
            convert_ms1 = False
        elif x == "--skip-ms2":
            print("<<<<: skipping ms2")
            convert_ms2 = False
        elif x[-1] == "*":
            for f in glob.glob(x):
                if is_valid_timstof_dir(f):
                    dirs_to_analyze.append(f)
        elif is_valid_timstof_dir(x):
            dirs_to_analyze.append(x)
        elif x.startswith("--output-path"):
            output_path = x[i+1]

    start_time = time.process_time()
    for timstof_path in dirs_to_analyze:
        place_high = 3
        t_path = timstof_path
        if t_path.endswith("/"):
            t_path = timstof_path[:len(timstof_path)-1]
        name = os.path.basename(t_path).split('.')[0]
        ms2_file_name = name + '_nopd.ms2'
        ms1_file_name = name + '_nopd.ms1'

        ms2_file_name = os.path.join(timstof_path, ms2_file_name)
        ms1_file_name = os.path.join(timstof_path, ms1_file_name)

        print(timstof_path)
        if convert_ms2:
            print(ms2_file_name)
        if convert_ms1:
            print(ms1_file_name)

        output = timstof_path + os.path.sep + ms2_file_name
        run_timstof_conversion(timstof_path, ms2_file_name)
    duration = start_time - time.process_time();
    print("duration is {}".format(duration))
