prm_ms2_scan_template = "S\t{frame:06d}\t{frame:06d}\t{monoisotopic_mz:.4f}\n" \
                        "I\tTIMSTOF_Prm_Frame_ID\t{frame}\n" \
                        "I\tTIMSTOF_Target_ID\t{target}\n" \
                        "I\tTIMSTOF_Scan_Begin\t{scan_begin}\n" \
                        "I\tIsolationMz\t{isolation_mz:.4f}\n" \
                        "I\tIsolationWidth\t{isolation_width:.4f}\n" \
                        "I\tCollisionEnergy\t{collision_energy:.4f}\n" \
                        "I\tOneOverK0\t{one_over_k0:.4f}\n" \
                        "I\tRetTime\t{ret_time:.4f}\n" \
                        "Z\t{cs:d}\t{m:.4f}\n"


class PrmTarget:
    def __init__(self, target_id: int, external_id: str, ret_time: float, one_over_k0: float, monoisotopic_mz: float,
                 charge: int, description: str):
        self.target_id = target_id
        self.external_id = external_id
        self.ret_time = ret_time
        self.one_over_k0 = one_over_k0
        self.monoisotopic_mz = monoisotopic_mz
        self.charge = charge
        self.description = description


class PrmFrame:
    def __init__(self, frame: int, scan_num_begin: int, scan_num_end: int, isolation_mz: float, isolation_width: float,
                 collision_energy: float, target: int):
        self.isolation_width = isolation_width
        self.target = target
        self.collision_energy = collision_energy
        self.isolation_mz = isolation_mz
        self.scan_num_end = scan_num_end
        self.scan_num_begin = scan_num_begin
        self.frame = frame


def create_prm_ms2_scan(prm_target: PrmTarget, prm_frame: PrmFrame):
    return prm_ms2_scan_template.format(frame=prm_frame.frame, monoisotopic_mz=prm_target.monoisotopic_mz,
                                        target=prm_target.target_id, scan_begin=prm_frame.scan_num_begin,
                                        one_over_k0=prm_target.one_over_k0, ret_time=prm_target.ret_time,
                                        collision_energy=prm_frame.collision_energy, isolation_mz=prm_frame.isolation_mz,
                                        isolation_width=prm_frame.isolation_width, cs=prm_target.charge,
                                        m=prm_frame.isolation_mz)