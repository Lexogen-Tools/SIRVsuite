
import matplotlib.pyplot as plt
        
y = [cnt_data[ctrl] for ctrl in sorted(self.conc.keys())]
cnt = float(sum(y))
y = [y_i/cnt if cnt>0 else 1e-100 for y_i in y]

x = [float(self.conc[ctrl][ercc_spike_in]) for ctrl in sorted(self.conc.keys())]
cnt = float(sum(x))
x = [x_i/cnt for x_i in x]
fig = plt.figure()
fig.gca().loglog(x,y, "rx")
x = np.array(x)
y = np.array(y)
ymin = 10**math.floor(math.log10(min(list(x[y>0])+list(y[y>0]))))
ymax = 10**math.ceil(math.log10(max(list(x[y>0])+list(y[y>0]))))
fig.gca().set_xlabel("Theoretical concentration")
fig.gca().set_ylabel("Measured concentration")
fig.gca().set_xlim(ymin, ymax)
fig.gca().set_ylim(ymin, ymax)
corr = scipy.stats.pearsonr([math.log10(x[i]) for i in range(len(x)) if x[i]!=0 and y[i]!=0],
                          [math.log10(y[i]) for i in range(len(y)) if x[i]!=0 and y[i]!=0])[0]
fig.gca().title.set_text("Sample %s (R=%.3f)"%
                         (helper.get_3tag_value(output[0]["tags"], "analysis", "alias"), corr))
plot_path = helper.pathjoin(output_dir, "ERCC.png")
fig.savefig(plot_path)
plt.close(fig)
data_path = os.path.join(output_dir, "spikein_stats.tsv")
with open(data_path,"w") as f:
    f.write("Alias\tERCC\n")
    f.write("%s\t%f\n"%(helper.get_3tag_value(output[0]["tags"], "analysis", "alias"), corr))

