
import curie as ci

def mca_to_spe_maestro_style(
    in_path,
    out_path=None,
    spec_id="No sample description was entered.",
    spec_rem_lines=None,
    meas_live=None,
    meas_real=None,
    date_str=None,
    indent=7,
    default_cal=(0.0, 1.0, 0.0)  # <- alltid skriv $MCA_CAL:, bruk denne hvis ingen finnes
):
    from datetime import datetime

    with open(in_path, "r", encoding="latin-1", errors="ignore") as f:
        txt = f.read()

    def _norm_lines(s):
        return [ln.strip() for ln in s.replace("\r\n","\n").replace("\r","\n").split("\n")]

    def _extract_block(joined, tag):
        start = joined.lower().find(f"<<{tag.lower()}>>")
        if start < 0:
            return None
        start += len(tag) + 4
        nxt = joined.find("<<", start)
        return joined[start:].strip() if nxt < 0 else joined[start:nxt].strip()

    def _parse_amptek_mca(s):
        lines = _norm_lines(s)
        joined = "\n".join(lines)
        data_block = _extract_block(joined, "DATA")
        if not data_block:
            return None
        dl = [l for l in _norm_lines(data_block) if l]
        counts = []
        if dl:
            try:
                n_decl = int(float(dl[0]))
                rest = [int(float(x)) for x in dl[1:]]
                counts = rest if n_decl == len(rest) else [int(float(x)) for x in dl]
            except:
                counts = [int(float(x)) for x in dl]

        pmca = _extract_block(joined, "PMCA SPECTRUM") or _extract_block(joined, "PMCA SPECTRUM DATA")
        ltime = rtime = None
        if pmca:
            for ln in _norm_lines(pmca):
                up, num, buf = ln.upper(), None, []
                for ch in ln:
                    if ch.isdigit() or ch in "+-.":
                        buf.append(ch)
                    elif buf:
                        try: num = float("".join(buf))
                        except: num = None
                        break
                if num is None and buf:
                    try: num = float("".join(buf))
                    except: num = None
                if num is None: continue
                if "LIVE" in up: ltime = num
                elif "REAL" in up: rtime = num

        # Amptek CALIBRATION
        calib = None
        cal = _extract_block(joined, "CALIBRATION")
        if cal:
            vals, tmp = [], ""
            for ch in cal:
                if ch.isdigit() or ch in "+-.": tmp += ch
                else:
                    if tmp:
                        try: vals.append(float(tmp))
                        except: pass
                        tmp = ""
            if tmp:
                try: vals.append(float(tmp))
                except: pass
            if len(vals) >= 3: calib = (vals[0], vals[1], vals[2])
            elif len(vals) == 2: calib = (vals[0], vals[1], 0.0)

        return {"counts": counts, "ltime": ltime, "rtime": rtime, "calib": calib}

    def _parse_plain(s):
        counts = []
        for ln in _norm_lines(s):
            if not ln: continue
            if all(ch.isdigit() or ch in "+-. " for ch in ln):
                try: counts.append(int(float(ln)))
                except: pass
        return counts if counts else None

    parsed = _parse_amptek_mca(txt)
    if parsed:
        counts, ltime, rtime, calib = parsed["counts"], parsed["ltime"], parsed["rtime"], parsed["calib"]
    else:
        counts, ltime, rtime, calib = _parse_plain(txt), None, None, None
        if not counts:
            raise ValueError("Ukjent .mca-format.")

    if meas_live is not None: ltime = float(meas_live)
    if meas_real is not None: rtime = float(meas_real)
    if date_str is None:
        date_str = datetime.now().strftime("%m/%d/%Y %H:%M:%S")  # US-format som Curie forventer

    if spec_rem_lines is None:
        spec_rem_lines = ["DET# 2","DETDESC# DETECTOR 2","AP# Maestro Version 7.01"]

    # Skriv alltid $MCA_CAL:
    if calib is None:
        calib = default_cal

    lines = []
    lines += ["$SPEC_ID:", spec_id]
    lines += ["$SPEC_REM:", *spec_rem_lines]
    lines += ["$DATE_MEA:", date_str]
    lines += ["$MEAS_TIM:", f"{(ltime or 0):.0f} {(rtime or 0):.0f}"]
    lines += ["$MCA_CAL:", f"{calib[0]:.8g} {calib[1]:.8g} {calib[2]:.8g}"]
    lines += ["$DATA:", f"0 {len(counts)-1}"]
    pad = " " * max(0, int(indent))
    lines += [f"{pad}{int(c)}" for c in counts]
    lines.append("")
    spe_text = "\n".join(lines)

    if out_path:
        with open(out_path, "w", encoding="ascii", errors="ignore") as f:
            f.write(spe_text)
        return out_path
    return spe_text



# new_spe_file = '/Users/hannah.ekeberg/Downloads/UiO_September2025-3/Calibration/CH10202025_Xray_Am241_p4_converted.Spe'
# old_mca_file = '/Users/hannah.ekeberg/Downloads/UiO_September2025-3/Calibration/CH10202025_Xray_Am241_p4.mca'
# mca_to_spe_maestro_style(old_mca_file, new_spe_file, date_str=None)
new_spe_file = '/Users/hannah.ekeberg/Downloads/UiO_September2025-3/Calibration/CE10182025_Xray_Co57_HC6928_p4_converted.Spe'
old_mca_file = '/Users/hannah.ekeberg/Downloads/UiO_September2025-3/Calibration/CE10182025_Xray_Co57_HC6928_p4.mca'
mca_to_spe_maestro_style(old_mca_file, new_spe_file, date_str=None)

# sp = ci.Spectrum(new_spe_file)
# sp.plot()