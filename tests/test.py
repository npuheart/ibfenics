import ibfenics

iu = ibfenics.UserIU(g=1,cm=1)



time_manager = ibfenics.io.TimeManager(1, 100, 20)


for i in range(1000):
    if time_manager.should_output(i):
        print(f"output {i}")