import numpy as np
import taichi as ti
import os
import imageio
import argparse
from tqdm import tqdm

class Viewer:
    def __init__(self, demo_name, resolution, frame_size):
        self.demo_name=demo_name
        self.resolution=resolution
        self.frame_size=frame_size
        self.frame_count=0

    def read_particles_file(self, filename):
        with open(filename, 'rb') as f:
            mm = np.memmap(f, dtype='<u4', mode='r')
            count = mm[0]
            positions = mm[1:1+count*2].view('<f4').reshape(-1, 2)
            radius = mm[1+count*2:1+count*2+1].view('<f4')[0]
            positions = positions * np.array([(self.resolution-4)/self.resolution, (self.resolution-4)/self.resolution/2]) + 0.5
        
        return count,positions, np.ceil(radius*self.frame_size*0.7)

    def draw(self):
        diction = "./build/windows/x64/release/output/"

        with open(diction + "frame_count.txt", "r") as f:
            self.frame_count = int(f.read())

        collider_particles_num, collider_particles_position, collider_particles_radius=self.read_particles_file(diction+f"results/0/collider.out")

        gui = ti.GUI(self.demo_name, (self.frame_size, self.frame_size*2))
        gui.background_color = 0xFFFFFF

        fluid_particles_position = []
        DEM_particles_position = []
        DEM_particles_radius = 0.
        fluid_particles_radius = 0.

        for i in tqdm(range(self.frame_count)):
            
            data_dict = diction + f"results/{i}/"
            _, temp_fluid_particles_position, fluid_particles_radius=self.read_particles_file(data_dict+"particles.out")
            _, temp_DEM_particles_position, DEM_particles_radius=self.read_particles_file(data_dict+"DEMparticles.out")
            fluid_particles_position.append(temp_fluid_particles_position)
            DEM_particles_position.append(temp_DEM_particles_position)
 
        for i in range(self.frame_count):
            gui.circles(fluid_particles_position[i], radius=fluid_particles_radius, color=0x0000FF)
            gui.circles(DEM_particles_position[i], radius=DEM_particles_radius, color=0xFFFF00)
            gui.circles(collider_particles_position, radius=collider_particles_radius, color=0x9F9F9F)
            filename = f'./results/{self.demo_name}/img/frame_{i}.png'
            gui.show(filename)


    def png_to_gif(self, input_dir, output_file, fps=30):
        filenames = sorted(
            [f for f in os.listdir(input_dir) if f.startswith("frame_") and f.endswith(".png")],
            key=lambda x: int(x.split("_")[1].split(".")[0])
        )[1:self.frame_count]
        
        with imageio.get_writer(output_file, mode='I', fps=fps, loop=0) as writer:
            for filename in filenames:
                image = imageio.imread(os.path.join(input_dir, filename))
                writer.append_data(image)
        
        print(f"GIF has been saved to {output_file}")

    def view(self):
        self.draw()
        self.png_to_gif(f"./results/{self.demo_name}/img", f"./results/{self.demo_name}/output.gif", fps=24)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='viewer')

    parser.add_argument('-t', '--demo_name', type=str,default='falling', help='demo_name')
    parser.add_argument('-s',  '--resolution',type=int, default=64, help='simulation resolution')
    parser.add_argument('-f', '--frame_size',type=int, default=400, help='frame size')

    args = parser.parse_args()

    frame_size = args.frame_size
    resolution=args.resolution
    demo_name =args.demo_name
    os.makedirs(f"./results/{demo_name}", exist_ok=True)
    os.makedirs(f"./results/{demo_name}/img", exist_ok=True)
    viewer = Viewer(demo_name,resolution,frame_size)
    viewer.view()