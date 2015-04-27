set table "shortest_path.g2.table"; set format "%.5f"
set samples 25; plot [x=1.0:7.0] (4.0/6.0)*x+1.0/3.0 + 0.05*(sin(1.5*x)**2+1)*(x-1.0)*(x-7.0)
