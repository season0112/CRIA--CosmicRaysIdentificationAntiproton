import matplotlib.pyplot as plt

def Plot_TimeDependentMeasuringTime(Time6B, MeasuringTime):

    fig_mine, ax = plt.subplots()

    #plt.bar(Time6B, MeasuringTime, width=datetime.fromtimestamp(6998400/2))
    plt.bar(Time6B, MeasuringTime, width=162, color='White', edgecolor='blue', linewidth=5)  #The unit for bar width on a date x axis is days. 6*27=162

    plt.hlines(13996800, Time6B[0]-(Time6B[1]-Time6B[0])/2, Time6B[-1]+(Time6B[1]-Time6B[0])/2, linestyles = "solid", colors='red', linewidth=7) # 27*6*24*60*60=13996800,

    #plt.vlines(Time6B[0]-(Time6B[1]-Time6B[0])/2 , 0, 13996800, colors = "black", linestyles = "dashed", linewidth=7)
    #plt.vlines(Time6B[-1]+(Time6B[1]-Time6B[0])/2, 0, 13996800, colors = "black", linestyles = "dashed", linewidth=7)

    plt.text(0.80, 0.85, "162 days", transform = ax.transAxes, fontsize=60, color='red')

    ax.axes.set_xlim([Time6B[0]-(Time6B[1]-Time6B[0])/2, Time6B[-1]+(Time6B[1]-Time6B[0])/2])
    ax.ticklabel_format(style='sci', scilimits=(0,0), axis='y', useMathText=True)
    ax.yaxis.get_offset_text().set_fontsize(60)

    plt.ylabel("T / (s)", horizontalalignment='right', y=1.0)

    plt.savefig( "TimeDependentMeasuringTime.pdf")
    plt.close()


