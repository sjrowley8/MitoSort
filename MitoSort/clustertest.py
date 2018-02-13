import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

with open('/home/m169420/cluster_test.txt', 'w') as f:
    f.write(str(plt.style.available))
    f.close()
