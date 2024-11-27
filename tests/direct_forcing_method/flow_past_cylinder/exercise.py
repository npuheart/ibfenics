#定义Dog类，类中的函数都称为方法，唯一不同是调用方法不同
# class Dog():
#     """ 一次模拟小狗的简单尝试"""
    
#     def __init__(self, name, age): 
#         """初始化属性 name和age """
#         self.name = name
#         self.age = age #以self为前缀的变量都可供类中的所有方法使用，我们可通过类的任何实例来访问这些变量，该变量可通过实例访问的变量称为属性
        
#     def sit(self):
#         """ 模拟小狗被命令时蹲下"""
#         print(self.name.title() + " is now sitting.")
        
#     def roll_over(self):
#         """模拟小狗被命令时打滚"""
#         print(self.name.title() + " rolled over!")
#     #Dog类还定义了另外两个方法，sit()和roll_over().这些方法不需要额外的信息，如名字和年龄，因此只需要一个形参。后面将创建的实例能够访问这些方法，换句话说，它们都会蹲下和
#     #打滚。
# #创建一个具体的实例
# dog = Dog("willie", 18) #python 使用实参'wille'和18调用Dog类中的方法_init_().该方法创建一个表示特定小狗的示例，并使用我们提供的值来设置属性name和age.
# #虽然_init_并未显式地包含return语句，但是python自动返回一个表示这条小狗的实例。将这个实例存储在变量my_dog中。 大写字母表示的是类，小写字母表示的是创建的实例。
# you_dog = Dog("lucy", 3 ) #即使给定第二条小狗指定同样的名字和年龄，Python依然会根据Dog类创建另一个实例
# #如何应用类中的方法
# dog.sit()
# dog.roll_over()  #使用句点表示法来调用Dog类中定义的任何方法
# dog.name   # 访问实例的属性，使用句点表示法这个就是访问类的属性
# print("My dog's name is " + dog.name.title()+ '.')
# print("Your dog is "+ you_dog.name.title()+ ".")
# #注意类的命名约定很有用：通常可认为首字母大写的名称（如Dog)指的是类，而小写的名称（如dog)指的是根据类创建的实例


# class Restaurant():
#     servant = 1
    
#     def __init__(self, restaurant_name, cuisine_type):
        
#         self.restaurant_name = restaurant_name
#         self.cuisine_type = cuisine_type
#         self.number_served = 0
     
#     def describe_restaurant(self):
#         print(self.restaurant_name.title() + " is call Pengfei restaurant")
#         print(self.cuisine_type.title())
        
#     def open_restaurant(self):
#         print("The restaurant is open for business as usual.")
    
#     def set_number_served(self, number):
#         self.number_served = number
        
#     def increment_number_served(self, number):
#         self.number_served += number
#         print(self.restaurant_name.title() + " has " + str(self.number_served) + ".")

# Pengfei = Restaurant("pengfei", "chuangcai")
# Pengfei.describe_restaurant()
# # print(Pengfei.number_served)
# Pengfei.number_served = 23
# print(Pengfei.number_served)
# Pengfei.set_number_served(23)
# print(Pengfei.number_served)
# Pengfei.increment_number_served(2)


# #类编好了，大部分时间都将花在使用根据类创建的实例上。

# class Car():
#     """ 一次模拟汽车的简单尝试 """
    
#     def __init__(self, make, model, year): #定义了方法_init_(),这个方法第一个形参为self,另外还包括三个形参 make, model和year
#         """初始化描述汽车的属性"""#方法_init_()接受这些形参的值，并将它们存储在根据这个类创建的实例的属性中。创建新的Car实例时，我们需要指定其制造商、型号、和生产年份。
#         self.make = make
#         self.model = model
#         self.year = year
    
#     def get_descriptive_name(self): #定义了一个铭文get_descriptive_name()的方法，使用属性year,make和model创建一个对汽车进行描述的字符串，无需打印每个属性的值
#         """返回整洁的描述信息"""  #为了在这个方法中访问属性的值，我们使用self.make、self.model 和self.year.
#         long_name = str(self.year) + ' ' +self.make + '  ' + self.model
#         return long_name.title()

# my_new_car = Car('audi', 'a4', 2016)  #根据Car类创建了一个实例，并将其存储到变量my_new_car中。
# print(my_new_car.get_descriptive_name()) 

# class Car():
#     def __init__(self, make, model, year):
#         """ 初始化描述汽车的属性 """  #类中的每个属性必须有初始值，哪怕这个值是0或空字符串。有些情况下，如设置默认值时，在方法_init_()内指定这种初始值是可行的；如果你对某个属性这样做了，无需包含为它提供初始值的形参。
#         self.make = make
#         self.model = model
#         self.year = year 
#         self.odometer_reading = 0  #这里添加一个名为odometer_reading 的属性，其初始值为0.
    
#     def update_odometer(self, mileage):
#         """将里程表读数设置为指定的值"""
#         self.odometer_reading = mileage   #对这个类唯一的修改是在Car类上添加了update_odometer().这个方法接受了一个里程值，并将其存储在self.odometer_reading中。
    
#     def get_descriptive_name(self):
#         """返回整洁的描述性信息"""
#         long_name = str(self.year) + ' '+ self.make + ' ' + self.model
#         return long_name.title()
    
#     def read_odometer(self):  #添加了一个铭文read_odometer()的方法，用来读取汽车的里程表
#         """打印一条指出汽车里程的消息"""
#         print("This car has " + str(self.odometer_reading) + " miles on it.")
        
#     def update_odometer(self, mileage):
#         """将里程表读数设置为指定值
#             禁止将里程表读数往回调
#         """
#         if  mileage >= self.odometer_reading: 
#             self.odometer_reading = mileage
#         else:
#             print("You can't roll back an odometer!")
            
#     def increment_odometer(self, miles): #新增的方法increment_odometer()接受一个单位为英里的数字，并将其加入到self.odometer_reading中。
#         """ 将里程表读数增加指定的量"""
#         self.odometer_reading += miles
        
# my_new_car = Car('audi', 'a4', 2016)
# print(my_new_car.get_descriptive_name())
# my_used_car = Car('subaru', 'outback', 2013)
# print(my_used_car.get_descriptive_name())
# my_new_car.odometer_reading = 23  #使用句点表示法来直接访问并设置汽车的属性odometer_reading.这行代码让Python在实例my_new_car中找到属性odometer_reading,并将该属性的值设置为23.
# my_new_car.update_odometer(23)
# my_new_car.read_odometer()  

#编写类时，并非总是要从空白开始，如果你要编写的类是一个现成类的特殊版本，可使用继承。一个类继承另一个类时，它将自动获得
#另一个类的所有属性和方法；原有的类称为父类，而新类称为子类。子类继承其父类的所有属性和方法，同时还可以定义自己的属性和方法。
#创建子类的实例时，Python首先需要完成的任务是给父类的所有属性赋值。为此，子类的方法_init_()需要父类。

#下面来模拟电动汽车。电动汽车是一种特殊的汽车，因此可以在前面创建的Car类的基础上创建新的类ElectricCar，这样只需为电动汽车特有的属性和行为编写代码

class Car(): #父类必须包含在当前文件中，且位于子类前面。
    """一次模拟汽车的简单尝试"""
    
    def __init__(self, make, model, year):
        self.make = make
        self.model = model 
        self.year = year
        self.odometer_reading = 0
        
    def get_descriptive_name(self):
        long_name = str(self.year) + ' ' + self.make + ' ' + self.model
        return long_name.title()
    
    def read_odometer(self):
        print("This car has " + str(self.odometer_reading) + " miles on it.")
    
    def update_odometer(self, mileage):
        if mileage >= self.odometer_reading:
            self.odometer_reading = mileage
        else:
            print("You can't roll back an odometer!")
    
    def increment_odometer(self, miles):
        self.odometer_reading += miles

class ElectricCar(Car): #定义了子类ElectricCar。定义子类时，必须在括号中指定父类的名称。_init_()接受创建Car实例所需的信息。
    """电动汽车的独特之处"""
    def __init__(self, make, model, year):
        """初始化父类的属性，再初始化电动车特有的属性"""
        super().__init__(make, model, year) #super()是一个特殊函数，帮助python将父类和子类关联起来。这里让python调用ElectricCar的父类方法_init_(),让
#ElectricCar 实例包含父类的所有属性。父类也称为超类（superclass),名称super因此而得名。
        self.battery_size = 70  #这里加入了新的属性，并设置了初始值
    
    def describe_battery(self): #这里加入了describe_battery的新方法，它打印有关电瓶的信息。
        """打印一条描述电瓶容量的消息"""
        print("This car has a " + str(self.battery_size) + "-kwh battery.")
        
    def fill_gas_tank():
        """电动汽车没有邮箱"""
        print("This car doesn't need a gas tank!") #如果调用电动车中的fill_gas_tank(),python将忽略Car类中的方法fill_gas_tank(),转而运行上述的代码。继承精华，弃其糟粕。
        
my_tesla = ElectricCar('tesla', 'models', 2016) #创建ElectricCar类的一个实例，并将其存储在变量my_tesla中。
print(my_tesla.get_descriptive_name())
my_tesla.describe_battery()
#对于父类的方法，只要它不符合子类模拟的实物行为，可对其重写。可在子类中定义这样一个方法，即它与要重写的父类方法同名。这样，python 将不会考虑这个父类方法，而
#只灌注在子类中定义的相应的方法。

#使用代码模拟实物时，可能会发现自己给类添加的细节越来越多；属性和方法清单以及文件都会越来越长。

