import timsdata, sqlite3, sys, time, os
import numpy as np, matplotlib.pyplot as plt
import sqlite3
from sqlite3 import Error

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


if __name__ == '__main__':
    if len(sys.argv)==1:
        print "Usage: extract_msn_nopd [source data directory (.d)] [target directory for output]"
    else:
        analysis_dir = sys.argv[1]

    start_time = time.clock()

    if sys.version_info.major == 2:
        analysis_dir = unicode(analysis_dir)
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

    precursor_array = np.array(precursor_list)  # 'ID', 'LargestPeakMz', 'AverageMz', 'MonoisotopicMz', 'Charge', 'ScanNumber', 'Intensity', 'Parent'
    frame_parent_dict = msms_frame_parent_dict(all_frame)

    # precursor_df = pd.DataFrame(data=precursor_array, index=precursor_array[:, 0], columns=['ID', 'LargestPeakMz', 'AverageMz', 'MonoisotopicMz', 'Charge', 'ScanNumber', 'Intensity', 'Parent'])

    parent_frame_array = np.array(precursor_array[:, 7])

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
    os.chdir(sys.argv[2])
    with open(ms2_file_name, 'wb') as output_file:
        output_file.write(ms2_header)
        progress = 0
        last_parent_frame = 0
        for each_msms_data in msms_data:
            frame_id, scan_begin, scan_end, precursor_mass, isolation_width, collision_energy, precursor_seq = each_msms_data
            # print frame_id,scan_begin, scan_end, precursor_mass
            frame_id_int = int(frame_id)
            scan_begin_int = int(scan_begin)
            scan_end_int = int(scan_end)
            scans = td.readScans(frame_id_int, scan_begin_int, scan_end_int)
            index_intensity = np.concatenate(scans, axis=1)
            mass_array = td.indexToMz(frame_id_int, index_intensity[0])
            mass_intensity = np.around(np.array(zip(mass_array, index_intensity[1])), decimals=4)
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
            for each_pair in sorted_mass_intensity:
                mass, intensity = each_pair
                if mass == last_mass:
                    total_intensity += intensity
                else:
                    output_file.write("%.4f %.1f\n" % (last_mass, enhance_signal(total_intensity)))
                    total_intensity = intensity
                last_mass = mass

            output_file.write("%.4f %.1f\n" % (last_mass, enhance_signal(total_intensity)))

            progress += 1
            if progress % 5000 == 0:
                print "progress %.1f%%" % (float(progress) / len(msms_data) * 100), time.clock() - start_time
