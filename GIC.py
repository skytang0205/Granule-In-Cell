import os
import argparse
from viewer import Viewer

parser = argparse.ArgumentParser(description='viewer')

parser.add_argument('-t', '--demo_name', type=str,default='falling', help='demo_name')
parser.add_argument('-s',  '--resolution',type=int, default=64, help='simulation resolution')
parser.add_argument('-f', '--frame_size',type=int, default=400, help='frame size')
parser.add_argument('-R', '--Radius',type=float, default=0.5, help='sand Radius')
parser.add_argument('-r', '--rate', type = int, default=50, help='Frame rate')
parser.add_argument('-e', '--end', type = int, default=200, help='End frame')

args = parser.parse_args()

frame_size=args.frame_size
resolution=args.resolution
demo_name =args.demo_name
Radius=args.Radius
rate=args.rate
end=args.end
os.makedirs(f"./results/{demo_name}", exist_ok=True)
os.makedirs(f"./results/{demo_name}/img", exist_ok=True)

os.system("xmake b demo")
os.system(f"xmake r demo -t {demo_name} -r {rate} -e {end} -s {resolution} -R {Radius}")

myViewer = Viewer(demo_name,resolution,frame_size)
myViewer.view()