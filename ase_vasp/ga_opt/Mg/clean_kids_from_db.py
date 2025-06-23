from ase.db import connect

db = connect("gadb.db")
lo = list(db.select("pairing=True"))[0].id
hi = list(db.select())[-1].id
print(f"Deleting rows {lo} through {hi}")
db.delete(ids=list(range(lo,hi+1)))
print("Done")
