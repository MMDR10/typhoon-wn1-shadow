import json, math

with open("wn1-shadow/wn1_history.json") as f:
    data = json.load(f)

storms = {}
for r in data["records"]:
    name = r["storm"].replace("TROPICAL DEPRESSION ", "")
    storms.setdefault(name, []).append(r)

for name in ["CHAN-HOM", "FIFTEEN"]:
    recs = sorted(storms[name], key=lambda x: x["cycle"])
    errs = []
    
    print(f"\n{name}:")
    for i in range(1, len(recs)):
        prev, curr = recs[i-1], recs[i]
        dlat = curr["center_lat"] - prev["center_lat"]
        dlon = curr["center_lon"] - prev["center_lon"]
        bearing = math.degrees(math.atan2(dlon, dlat)) % 360
        
        phi = curr["wn1_phi"]
        err = abs(phi - bearing)
        if err > 180: err = 360 - err
        
        errs.append(err)
        print(f"  {curr['cycle']}: phi={phi:5.1f} bear={bearing:5.1f} err={err:5.1f}")
    
    print(f"  ---")
    print(f"  Mean: {sum(errs)/len(errs):.1f}")
    print(f"  Median: {sorted(errs)[len(errs)//2]:.1f}")
    print(f"  Samples: {len(errs)}")
