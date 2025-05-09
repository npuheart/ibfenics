

import swanlab
swanlab.login(api_key="VBxEp1UBe2606KHDM9264", save=True)
                                                                                                
swanlab.init(
    project="AFSI",
    experiment_name="p1",
    description="This is a test experiment",
    config={'learning-rate': 0.003},
)


# 记录指标
for i in range(10):
    swanlab.log({"loss": i, "acc": i-1})