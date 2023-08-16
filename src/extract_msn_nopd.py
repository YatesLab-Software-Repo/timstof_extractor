from enum import Enum

import timsdata, sqlite3, sys, time, os
import numpy as np, matplotlib.pyplot as plt
from scipy import constants
import sqlite3
from sqlite3 import Error
import glob
import time
import math
from datetime import datetime

from dia_frame import DiaFrame
from ms_string_templates import header_ms2_template, header_ms1_template, \
    dda_ms2_scan_template
from src.prm_data import PrmTarget, PrmFrame, create_prm_ms2_scan

place_high = 3
precursor_counter = 0
convert_ms2 = True
convert_ms1 = True
rename = False
version = "0.3.0"
batch_size = 1000

class DataType(Enum):
    DIA = 0
    DDA = 1
    PRM = 2
# no batching 6 m 16 sec
# at batch size 100  time was 5 m 7 s
# at batch size 500  time was 4 m 34 s
# at batch size 1000 time was  4 m 16.427 s
# at batch size 1500 time was 4 m 22.28 s
# at batch size 2000 time was 4 m 36.07 s

def K0toCCS(K0, q, m_ion, m_gas, T):
    mu = m_ion * m_gas / (m_ion + m_gas)
    T0 = 273.15
    p0 = 1.0132e5  # 1 atm
    N0 = p0 / (constants.k * T0)
    return (3.0 / 16.0) * (1 / N0) * np.sqrt(2 * np.pi / (mu * constants.k * T)) * (q / K0)


def enhance_signal(intensity):
    return (intensity ** 1.414 + 100.232) ** 1.414


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


def select_all_prm_frame_msms_info(conn):
    cur = conn.cursor()
    rows = cur.execute(
        'select Frame, ScanNumBegin, ScanNumEnd, IsolationMz, IsolationWidth, CollisionEnergy, Target from '
        'PrmFrameMsMsInfo')
    prm_frame_list = [
        PrmFrame(frame=Frame, scan_num_begin=ScanNumBegin, scan_num_end=ScanNumEnd, isolation_mz=IsolationMz,
                 isolation_width=IsolationWidth, collision_energy=CollisionEnergy, target=Target)
        for Frame, ScanNumBegin, ScanNumEnd, IsolationMz, IsolationWidth, CollisionEnergy, Target in rows]
    return prm_frame_list


def select_all_prm_targets(conn):
    cur = conn.cursor()
    rows = cur.execute('select Id, ExternalId, Time, OneOverK0, MonoisotopicMz, Charge, Description from PrmTargets')
    prm_target_list = {Id:
                           PrmTarget(target_id=Id, external_id=ExternalId, ret_time=Time, one_over_k0=OneOverK0,
                                     monoisotopic_mz=MonoisotopicMz,
                                     charge=Charge, description=Description)
                       for Id, ExternalId, Time, OneOverK0, MonoisotopicMz, Charge, Description in rows}

    return prm_target_list


def select_all_dia_frames(conn):
    cur = conn.cursor()
    rows = cur.execute('select DiaFrameMsMsWindows.windowgroup, frame, scannumbegin, scannumend, IsolationMz, '
                       'IsolationWidth, CollisionEnergy, Frames.time  from DiaFrameMsMsWindows join DiaFrameMsMsInfo on '
                       'DiaFrameMsMsWindows.windowgroup = DiaFrameMsMsInfo.windowgroup '
                       'join Frames on DiaFrameMsMsInfo.Frame = Frames.id')

    dia_list = [
        DiaFrame(window_group=window_group, frame=frame, scan_num_begin=scan_num_begin, scan_num_end=scan_num_end,
                 isolation_mz=isolation_mz, isolation_width=isolation_width, collision_energy=collision_energy,
                 time=time)
        for window_group, frame, scan_num_begin, scan_num_end, isolation_mz, isolation_width, collision_energy,
            time in rows]

    return dia_list


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
            for count, rows in enumerate(precursor_map[frame_id], prev_scan + 1):
                prec_id = int(rows[0])
                ms2_map[frame_id][prec_id] = count
            prev_scan += len(precursor_map[frame_id])
    return frame_id_ms1_scan_map, ms2_map


def build_dia_frame_id_ms1_scan_map(all_ms1_frame, frame_dia_map):
    prev_scan = 0
    frame_id_ms1_scan_map = {}
    ms2_map = {}
    for row in all_ms1_frame:
        frame_id = int(row[0])
        prev_scan += 1
        frame_id_ms1_scan_map[frame_id] = prev_scan
        if frame_id in frame_dia_map:
            for dia in frame_dia_map[frame_id]:
                prev_scan += 1
                ms2_map.setdefault(dia.window_group, {})[dia.scan_num_begin] = prev_scan
    return frame_id_ms1_scan_map, ms2_map


def create_mz_int_spectra(mz_int_arr):
    mz_arr = mz_int_arr[0]
    str_list = []
    int_arr = mz_int_arr[1]
    for j in range(0, len(mz_arr)):
        str_list.append("%.4f %.1f \n" % (mz_arr[j], int_arr[j]))
    return ''.join([row for row in str_list])


def create_dia_ms2_file(td, dia_list, ms2_file_name):
    with open(ms2_file_name, 'w') as output_file:
        progress = 0
        filler_num = 1
        for dia in dia_list:
            scan_num_arr = [scan_num for scan_num in range(dia.scan_num_begin, dia.scan_num_end)]
            one_over_k0_arr = td.scanNumToOneOverK0(dia.frame, scan_num_arr)
            index_intensity_arr = td.readScans(dia.frame, dia.scan_num_begin, dia.scan_num_end)
            index_intensity_carr = np.concatenate(index_intensity_arr, axis=1)
            index_arr = index_intensity_carr[0]
            mass_array = td.indexToMz(dia.frame, index_arr)
            index_mass_dict = dict(zip(index_arr, mass_array))
            spectra = []
            for scan_num_index, row in enumerate(index_intensity_arr):
                ook0 = one_over_k0_arr[scan_num_index]
                for ii in range(len(row[0])):
                    mz = index_mass_dict[row[0][ii]]
                    intensity = row[1][ii]
                    spectra.append((mz, intensity, ook0))
            spectra.sort(key=lambda tup: tup[0])
            scan_head = dia.get_scan_head(filler_num)
            output_file.write(scan_head)
            for mz, intensity, ook0 in spectra:
                output_file.write("{:.4f} {:.4f} {:.4f}\n".format(mz, intensity, ook0))
            filler_num += 1


def create_dda_ms2_file(td, precursor_list, ms2_scan_map, all_frame, date_now, ms2_file_name):
    with open(ms2_file_name, mode='wt', buffering=10) as output_file:
        first_scan = ms2_scan_map[int(precursor_list[0][7])][int(precursor_list[0][0])]
        last_scan = ms2_scan_map[int(precursor_list[-1][7])][int(precursor_list[-1][0])]
        ms2_header = header_ms2_template.format(version=version, date_of_creation=date_now, first_scan=first_scan,
                                                last_scan=last_scan)
        output_file.write(ms2_header)
        progress = 0
        row_retrieve_list = []
        prc_id_list = []
        parent_scan_dict = {}
        for row in precursor_list:
            prc_id, largest_preak_mz, average_mz, monoisotopic_mz, cs, scan_number, intensity, parent = row
            prc_id_int = int(prc_id)
            if monoisotopic_mz is not None and cs is not None and cs > 1:
                row_retrieve_list.append(row)
                prc_id_list.append(prc_id_int)

        row_batch = [row_retrieve_list[x:x + batch_size] for x in range(0, len(row_retrieve_list), batch_size)]
        prc_id_batch_list = [prc_id_list[x:x + batch_size] for x in range(0, len(prc_id_list), batch_size)]

        for i in range(0, len(row_batch)):
            rb = row_batch[i]
            prc_batch = prc_id_batch_list[i]
            mz_int_arr = td.readPasefMsMs(prc_batch)
            for row in rb:
                prc_id, largest_preak_mz, average_mz, monoisotopic_mz, cs, scan_number, intensity, parent = row
                prc_id_int = int(prc_id)
                spectra_arr = mz_int_arr[prc_id_int]
                if len(spectra_arr[0]) > 0:
                    prc_mass_mz = float(monoisotopic_mz)
                    prc_mass = (prc_mass_mz * cs) - (cs - 1) * 1.007276466
                    parent_index = int(parent)
                    scan_id = ms2_scan_map[parent_index][prc_id_int]
                    rt_time = float(all_frame[parent_index - 1][1])
                    k0 = td.scanNumToOneOverK0(parent_index, [scan_number])[0]
                    #  k0 = parent_scan_k0_matrix[parent_index][scan_number]

                    scan_head = dda_ms2_scan_template.format(scan_id=scan_id, prc_mass_mz=prc_mass_mz,
                                                             parent_index=parent_index, prc_id=prc_id_int,
                                                             ret_time=rt_time, k0=k0, cs=int(cs), prc_mass=prc_mass)
                    output_file.write(scan_head)
                    mz_arr = spectra_arr[0]
                    int_arr = spectra_arr[1]
                    for j in range(0, len(mz_arr)):
                        output_file.write("%.4f %.1f \n" % (mz_arr[j], int_arr[j]))

                progress += 1
                if progress % 5000 == 0:
                    print("progress ms2: {:.1f}% {:.1f}".format((float(progress) / len(precursor_list) * 100),
                                                                time.time() - start_time))


def create_prm_ms2_file(td, ms2_file_name):
    with open(ms2_file_name, "w") as ms2_file:
        prm_target_map = select_all_prm_targets(td.conn)
        prm_frame_list = select_all_prm_frame_msms_info(td.conn)
        prm_frame_list.sort(key=lambda x: x.frame)
        for prm_frame in prm_frame_list:
            prm_target = prm_target_map[prm_frame.target]
            index_intensity_arr = td.readScans(prm_frame.frame, prm_frame.scan_num_begin, prm_frame.scan_num_end)
            if len(index_intensity_arr) > 0:
                scan_num_arr = [scan_num for scan_num in range(prm_frame.scan_num_begin, prm_frame.scan_num_end)]
                one_over_k0_arr = td.scanNumToOneOverK0(prm_frame.frame, scan_num_arr)
                index_intensity_carr = np.concatenate(index_intensity_arr, axis=1)
                index_arr = list(set(index_intensity_carr[0]))
                mass_array = td.indexToMz(prm_frame.frame, index_arr)
                index_mass_dict = dict(zip(index_arr, mass_array))
                idx_intensity_map={}
                for scan_num_index, row in enumerate(index_intensity_arr):
                    ook0 = one_over_k0_arr[scan_num_index]
                    for ii in range(len(row[0])):
                        idx = row[0][ii]
                        mz = index_mass_dict[row[0][ii]]
                        intensity = row[1][ii]
                        entry = idx_intensity_map.get(idx)
                        if entry is None:
                            entry = [mz, intensity, ook0]
                        else:
                            entry[1] = entry[1] + intensity
                        idx_intensity_map[idx] = entry
                spectra = [entry for entry in idx_intensity_map.values()]
                spectra.sort(key=lambda tup : tup[0])
                ms2_scan_head = create_prm_ms2_scan(prm_target, prm_frame)
                ms2_file.write(ms2_scan_head)
                for row in spectra:
                    ms2_file.write("{:.4f} {:.4f} {:.4f}\n".format(row[0], row[1], row[2]))


def run_timstof_conversion(input, ms2_file_path, ms1_file_name="", data_type=DataType.DDA):
    global place_high
    global precursor_counter
    analysis_dir = input

    td = timsdata.TimsData(analysis_dir)
    conn = td.conn

    d = str(datetime.now().strftime("%B %d, %Y %H:%M"))
    if data_type == DataType.PRM:
        create_prm_ms2_file(td, ms2_file_name=ms2_file_path)
    elif data_type == DataType.DIA:
        with conn:
            all_frame = select_all_Frames(conn)
            all_ms1_frames = [a for a in all_frame if a[4] == '0']
            dia_frame_list = select_all_dia_frames(conn)
            frame_dia_map = {}
            for dia in dia_frame_list:
                frame_dia_map.setdefault(dia.frame, []).append(dia)
            # build_dia_frame_id_ms1_scan_map(all_ms1_frames, frame_dia_map)
            create_dia_ms2_file(td, dia_frame_list, ms2_file_name=ms2_file_path)
    else:
        precursor_map = {}
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
        if convert_ms2:
            create_dda_ms2_file(td=td, precursor_list=precursor_list, all_frame=all_frame, date_now=d,
                                ms2_scan_map=ms2_scan_map, ms2_file_name=ms2_file_path)
        if convert_ms1:
            ms1_header = header_ms1_template.format(version)
            with open(ms1_file_name, 'w') as output_file:
                output_file.write(ms1_header)
                progress = 0
                prev_id = 0
                # scan_set = set()
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
                    if len(sorted_mass_intensity) > 0:
                        rt_time = 0 if i == 0 else all_ms1_frames[i - 1][1]
                        lines.append("S\t%06d\t%06d\n" % (scan_num, scan_num))
                        lines.append("I\tTIMSTOF_Frame_id\t{}\n".format(id))
                        lines.append("I\tRetTime\t%.2f\n" % float(rt_time))
                        for row in sorted_mass_intensity:
                            x_str = "%.4f %.1f %.4f \n" % (row[0], row[1], row[-2])
                            lines.append(x_str)
                    if len(lines) > 500_000:
                        output_file.writelines(lines)
                        lines = []

                    progress += 1
                    if progress % 5000 == 0:
                        print("progress ms1 %.1f%%" % (float(progress) / len(all_ms1_frames) * 100), time.process_time()
                              - start_time)
                output_file.writelines(lines)
                lines = []
    conn.close()
    if rename:
        for file in os.listdir(analysis_dir):
            if file == "analysis.tdf":
                tdf_new_name = ms2_file_path.replace(".ms2", ".tdf")
                os.rename(os.path.join(analysis_dir, file), tdf_new_name)
            if file == "analysis.tdf_bin":
                tdf_bin_new_name = ms2_file_path.replace(".ms2", ".tdf_bin")
                os.rename(os.path.join(analysis_dir, file), tdf_bin_new_name)


def is_valid_timstof_dir(path):
    if os.path.isdir(path):
        tdf_list = [f for f in os.listdir(path) if f.endswith(".tdf")]
        bin_list = [f for f in os.listdir(path) if f.endswith(".tdf_bin")]
        if len(tdf_list) == 1 and len(bin_list) == 1:
            return True
    return False


if __name__ == '__main__':
    dirs_to_analyze = set()
    print("Running timSTOFExtractor {}".format(version))
    if len(sys.argv) == 1:
        print("Usage: extract_msn_nopd [source data directory (.d)]")
        print("--skip-ms1: skips creating ms1 files")
        print("--skip-ms2: skips creating ms2 files")
    else:
        analysis_dir = sys.argv[1]
    idx_to_skip_set = set()
    data_type = DataType.DDA
    output_path = ""
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
                    dirs_to_analyze.add(f)
        elif x.startswith("-output-path"):
            output_path = sys.argv[i + 1]
            idx_to_skip_set.add(i + 1)
        elif x.startswith("--dia"):
            data_type = DataType.DIA
        elif x.startswith("--prm"):
            data_type = DataType.PRM

    for i, x in enumerate(sys.argv):
        if i not in idx_to_skip_set:
            if is_valid_timstof_dir(x):
                dirs_to_analyze.add(x)

    start_time = time.time()
    for timstof_path in dirs_to_analyze:
        place_high = 3
        t_path = timstof_path
        if t_path.endswith(os.path.sep):
            t_path = timstof_path[:len(timstof_path) - 1]
        name = os.path.basename(t_path).split('.')[0]
        ms2_file_name = name + '.ms2'

        print(timstof_path)
        if convert_ms2:
            print(ms2_file_name)
        if len(output_path) > 0:
            ms2_path = os.path.join(output_path, ms2_file_name)
        else:
            ms2_path = os.path.join(timstof_path, ms2_file_name)
        run_timstof_conversion(timstof_path, ms2_file_path=ms2_path, data_type=data_type)
    duration = (time.time() - start_time)
    min_dur = int(duration / 60)
    sec_dur = duration % 60
    print("duration is {:d} m {:.3f} s".format(min_dur, sec_dur))
