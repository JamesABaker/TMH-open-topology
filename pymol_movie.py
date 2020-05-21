import datetime

for n, i in enumerate(range(72)):
    cmd.turn( "y", 5.0 )
    cmd.ray(1000,1000)
    cmd.png(f"pymol_gifs/{datetime.datetime.now()}.png")
    print(f"Completed frame {n}")
