# INS_errors_model
A project, which include model errors of INS (inertial navigation sysytem). This model need to imitate operation of navigation devices (gyros, acsselerometrs). At the etrance arrives file with aircraft's position (coordinats in  navigation coordinate system North-Up-East) and velocities (.txt file). At the end model get file .txt with parametrs (velocities, coordinates), including errors appliances.

To use my prog need to add all py files in project, using Python 2.7 and run main.py. 

To run project, need to put file with telemetries aircrafts (simple file is listed) with name tt_from_ipo.txt. then run the prog and prog draw some grafics.

This project was created with algorithms, which was description in:
a. D. Titterton и J. Weston, «Strapdown Inertial Navigation Technology,» The Institution of Electrical Engineers, Reston, 2004.
b. Г. И. Емельянцев и А. П. Степанов, «Интегрированные инерциально-спутниковые системы ориентации и навигации,» ГНЦ РФ АО "Концерн "ЦНИИ "Электроприбор", СПб, 2016. (for eanglish need to translate)
c. А. А. Голован, Н. А. Парусников, Н. Б. Вавилова, М. Ю. Попеленский, О. Н. Богданов и А. А. Панев, «Математические модели уравнения ошибок БИНС и корректирующих измерений,» МГУ, Москва, 2009. (for english need to translate)

Project structire description:

main.py - header file which have basic algorothm
Bins_function.py - file which have all models of errors INS
plotting.py - file which have functions to plot the graphics 
read_from_file.py - file need to read txt files with telem
