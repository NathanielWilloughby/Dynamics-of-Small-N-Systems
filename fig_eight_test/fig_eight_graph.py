import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

x1, x2, x3 = np.loadtxt('xpositions.txt', unpack=True)
y1, y2, y3 = np.loadtxt('ypositions.txt', unpack=True)
Eerr, t = np.loadtxt('energy_error.txt', unpack=True)


fig, (ax_1, ax_2, ax_3) = plt.subplots(nrows=1, ncols=3, layout='tight', figsize=(18,3))
ax_1.plot(x1, y1)
ax_2.plot(x2, y2)
ax_3.plot(x3, y3)
ax_1.set_xlabel('Position along x-axis')
ax_2.set_xlabel('Position along x-axis')
ax_3.set_xlabel('Position along x-axis')
ax_1.set_ylabel('Position along y-axis')
ax_2.set_ylabel('Position along y-axis')
ax_3.set_ylabel('Position along y-axis')
ax_1.set_title('Particle 1')
ax_2.set_title('Particle 2')
ax_3.set_title('Particle 3')
plt.savefig(fname='fig_eight_trajectories.png')
plt.show()
plt.plot(t, Eerr)
plt.xlabel('Time')
plt.ylabel('Fractional Energy Error')
plt.savefig(fname='fig_eight_energy_error.png')
plt.show()

nframes = np.shape(t)[0]
animfig, ax = plt.subplots()
frames = []

for i in range(nframes):
    frame = ax.scatter((x1[i], x2[i], x3[i]), (y1[i], y2[i], y3[i]), color='blue')
    frames.append([frame])

ax.set_xlabel('Position along x-axis')
ax.set_ylabel('Position along y-axis')
animfig.suptitle('Figure of Eight')
anim = animation.ArtistAnimation(fig=animfig, artists=frames, repeat=True, interval=1)
anim.save(filename='fig_eight_animation.gif', writer='pillow')
plt.show()

    
    


