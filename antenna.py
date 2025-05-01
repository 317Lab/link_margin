from PyNEC import *
import numpy as np
class Wire:
    def __init__(self, x1,y1,z1,x2,y2,z2,rad,num_segs,tag_num):
        self.x1=x1
        self.y1=y1
        self.z1=z1
        self.x2=x2
        self.y2=y2
        self.z2=z2
        self.rad=rad
        self.num_segs=num_segs
        self.tag_num=tag_num
wire_list = []
wire_list.append(Wire(0.0,-0.413,0.0,0.0,0.0,0.0,0.002,19,1))
wire_list.append(Wire(0.0,0.0,0.0,0.0,0.413,0.0,0.002,19,2))
wire_list.append(Wire(0.0,-0.448,-0.294,0.0,0.0,-0.294,0.002,19,3))
wire_list.append(Wire(0.0,0.0,-0.294,0.0,0.448,-0.294,0.002,19,4))
wire_list.append(Wire(-0.413,0.0,0.0,0.0,0.0,0.0,0.002,19,5))
wire_list.append(Wire(0.0,0.0,0.0,0.413,0.0,0.0,0.002,19,6))
wire_list.append(Wire(-0.448,0.0,-0.294,0.0,0.0,-0.294,0.002,19,7))
wire_list.append(Wire(0.0,0.0,-0.294,0.448,0.0,-0.294,0.002,19,8))


context = nec_context()
geo = context.get_geometry()
for i in wire_list:
    geo.wire(tag_id=i.tag_num,segment_count=i.num_segs, xw1=i.x1,yw1=i.y1,zw1=i.z1,xw2=i.x2,yw2=i.y2,zw2=i.z2,rad=i.rad,rdel=1.0,rrad=1.0)
#geo.generate_cylindrical_structure(4,18)
context.geometry_complete(0)
print(type(geo))
#geo.wire(tag_id=1, segment_count = nr_segments, xw1=-0.5,yw1=0,zw1=0.5,xw2=0.5,yw2=0,zw2=0.5)
# define top wire 2
