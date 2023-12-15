import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import imageio

sns.color_palette("viridis", as_cmap=True)

#for i in range(1, 999):
#	#print(i, f'result{i}.dat')
#	df = pd.read_csv(f'result{i}.dat', sep=r'\s{2,}', header=None,  engine='python')
#	fig = sns.heatmap(df, xticklabels=False, yticklabels=False, cmap="Spectral")
#	plt.savefig(f'saved{i}.png')
#	plt.clf()
	
df = pd.read_csv(f'result999.dat', sep=r'\s{2,}', header=None,  engine='python')
fig = sns.heatmap(df, xticklabels=False, yticklabels=False, cmap="Spectral")
plt.savefig(f'saved999.png')
plt.clf()

df_1 = pd.read_csv(f'result1.dat', sep=r'\s{2,}', header=None,  engine='python')
df_500 = pd.read_csv(f'result500.dat', sep=r'\s{2,}', header=None,  engine='python')
df_999 = pd.read_csv(f'result999.dat', sep=r'\s{2,}', header=None,  engine='python')

fig, axs = plt.subplots(nrows = 3, ncols = 1)
fig.suptitle('Разрезы вдоль осей')

x = [i for i in range(len(df_1.iloc[50]))]
y = [i for i in range(len(df_1.iloc[:, 50]))]

axs[0].plot(x, df_1.iloc[50])
axs[0].plot(y, df_1.iloc[:, 50])
axs[0].set_ylabel('1')
axs[1].plot(x, df_500.iloc[50])
axs[1].plot(y, df_500.iloc[:, 50])
axs[1].set_ylabel('500')
axs[2].plot(x, df_999.iloc[50])
axs[2].plot(y, df_999.iloc[:, 50])
axs[2].set_ylabel('999')

plt.savefig('', dpi=300)

images = []

#for i in range(1, 999):
#	images.append(f'saved{i}.png')

#imageio.mimsave('animation.gif', [imageio.imread(image) for image in images], duration=0.1)

