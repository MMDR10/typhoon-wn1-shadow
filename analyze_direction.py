import json, math

with open("wn1-shadow/wn1_history.json") as f:
    data = json.load(f)

# Group by storm
storms = {}
for r in data["records"]:
    name = r["storm"].replace("TROPICAL DEPRESSION ", "")
    storms.setdefault(name, []).append(r)

for name, recs in storms.items():
    recs.sort(key=lambda x: x["cycle"])
    print(f"\n{'='*70}")
    print(f"  {name} ({len(recs)} samples)")
    print(f"{'='*70}")
    
    prev_bearing = None
    for i, r in enumerate(recs):
        cycle = r["cycle"]
        lat, lon = r["center_lat"], r["center_lon"]
        phi = r["wn1_phi"]
        ell = r["ellipt"]
        
        if i == 0:
            print(f"  {cycle} | {lat:5.1f}N {lon:5.1f}E | phi={phi:6.1f} | (first)    |        |          |          | ell={ell:.3f}")
            continue
            
        prev = recs[i-1]
        dlat = lat - prev["center_lat"]
        dlon = lon - prev["center_lon"]
        bearing = math.degrees(math.atan2(dlon, dlat)) % 360
        dist_km = math.sqrt((dlat*111)**2 + (dlon*111*math.cos(math.radians(lat)))**2)
        
        err = abs(phi - bearing)
        if err > 180: err = 360 - err
        
        dphi = phi - prev["wn1_phi"]
        if dphi > 180: dphi -= 360
        if dphi < -180: dphi += 360
        
        if prev_bearing is not None:
            dbearing = bearing - prev_bearing
            if dbearing > 180: dbearing -= 360
            if dbearing < -180: dbearing += 360
            db_str = f"{dbearing:+6.1f}"
        else:
            db_str = "     -"
        
        print(f"  {cycle} | {lat:5.1f}N {lon:5.1f}E | phi={phi:6.1f} | bear={bearing:5.1f} | err={err:5.1f} | d_phi={dphi:+6.1f} | d_bear={db_str} | ell={ell:.3f} | {dist_km:.0f}km")
        prev_bearing = bearing
