import ibfenics1

iu = ibfenics1.UserIU(g=1, cm=1)


time_manager = ibfenics1.io.TimeManager(1, 100, 20)


for i in range(1000):
    if time_manager.should_output(i):
        print(f"output {i}")
