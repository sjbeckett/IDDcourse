import marimo

__generated_with = "0.10.12"
app = marimo.App()


@app.cell(hide_code=True)
def _():
    #import some packages with functions to help run the notebook
    import marimo as mo #notebook package!
    import matplotlib.pyplot as plt   #helps with plotting
    import numpy as np    #helps with arrays
    from scipy.integrate import odeint   #helps with numerical integration of dynamical systems
    return mo, np, odeint, plt


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        # Simple intervention model
        #### Dr. Stephen Beckett, University of Maryland (beckett@umd.edu)
                **Welcome** to the simulation exercise section utilized in slides for Day 5 of BSCI439C/BIOL708F: "Infectious disease dynamics: a systems approach".

                Recall, the key point in working through this interactive notebook is not in learning how to code, which is challenging and goes beyond what we offer in this course, but to see how concepts we have discussed can be translated into models and simulations -- and the types of characteristics and dynamical outcomes these models describe.

        ## About this notebook
                This notebook was designed and written using a marimo.io notebook which is coded in python - a general computing programming language. The bonus of marimo is that it can be used as an interactive notebook environment with reactivity! This means that we can use objects such as sliders to change input parameters, that will automatically update outputs - such as figures showing simulation data. This makes the notebook easy to use -- even without needing to learn computer programming! **Remember**, you can download images by right clicking on them, or finding the 'Export output as PNG' in the ... of the appropriate figure cell.

        ## About this modeling exercise

                Today, we will use a model which is a slight modification of those we have seen previously - the SEIRD model. This model includes two infection stages, those who are exposed to a disease $E$ but not yet infectious, and those who are infectious $I$. Here, individuals are infectious for an average duration (1/$\gamma$) with a fraction $f$ dying and a fraction $(1-f)$ recovering from the infection, such that $f$ reflects the **infection fatality ratio** - the number of fatalities per infection. Today we will use a model formulation that scales such that our state variables are in units of people, rather than fractions of the population. The model is as follows:

                $$\begin{align}\frac{dS}{dt} =& -\beta \frac{SI}{N} \nonumber\\
                \frac{dE}{dt} =& \beta \frac{SI}{N} - \nu E \nonumber\\
                \frac{dI}{dt} =& \nu E - \gamma I \nonumber\\
                \frac{dR}{dt} =& (1-f)\gamma I\nonumber\\
                \frac{dD}{dt} =& f\gamma I.\nonumber\end{align}$$

        We will assume a population of $N=100,000$ people, one of whom is initially exposed to the infectious disease.

        Additionally, we will consider what the transmission rate, $\beta$ reflects: it is both about how often there may be close contact with an infectious individual, but also about how likely transmission is to occur given contact. Let $c$ be the average number of contacts that an individual is likely to have in a day; and $p$ be the chance that a contact could lead to infection. We can write this as:

        $$\beta = c.p$$

        We will use the above model to simulate an outbreak, and some potential responses. Remember, what follows contains many simplified assumptions regarding how a real outbreak, and responses to it, might unfold.

        ##Outbreak!

        The disease will spread unabated in the population until it is noticed. We will denote this time **$t_a = 30$ days** after the intial exposure. After this time there is potential for the population, public health bodies and goverment to react -- which has the potential to change the underlying properties of how the disease is spreading (rate parameters in our model), and ultimately the disease dynamics. Yet, we note not all available information up to day $t_a$ may be available to the population. In the event of an outbreak of unknown origin, deaths may serve as the best indicator, in the absence of testing. Additionally, can only track backwards to the time of first death -- and potentially to origin of that individuals symptoms; but this may not match the time the population was originally exposed to the disease. See plot below. Note, in real situations the number of attributed deaths could also be undercounted, especially if symptoms are compatible with other causes.

        For sake of simplicity, we will specify a changepoint in the model when parameter values change from one value to another -- the timing of this changepoint, $t_b$, is used to represent when the population reacts after they know there is a disease outbreak - this could be sometime after the outbreak alert $t_a$. Different types of reaction may have different consequences.
        """
    )
    return


@app.cell(hide_code=True)
def _(D00, N, np, plt, t00):
    fig = plt.figure(figsize=(8,4))
    plt.scatter(t00[0:28],np.floor(D00[0:28]),label='D obs')
    plt.vlines(30,ymin=-0.5,ymax=2*N,linestyles='dashed',colors='black') #represents today
    plt.xlabel('Time (days)')
    plt.ylabel('Cumulative deaths')
    plt.xlim([16.5,30.5])
    plt.ylim([-1,max(D00[0:28])+5])
    plt.legend(loc="upper left")
    plt.title(f"Observed cumulative death in outbreak of unknown origin")
    plt.gca()
    return (fig,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## Intervention:
        After this time there is potential for the population, public health bodies and goverment to react -- which has the potential to change the underlying properties of how the disease is spreading (rate parameters in our model), and ultimately the disease dynamics. For sake of simplicity, we will specify a changepoint in the model when parameter values change from one value to another -- the timing of this changepoint, $t_b$, is used to represent when the population reacts after they know there is a disease outbreak - this could be sometime after the outbreak alert $t_a$. Different types of reaction may have different consequences.

        ### 1) Take no action:
        In the absence of any reaction -- with no behavior change, no masking or changing of plans, no change even in response to news of the disease, then the underlying parameters will stay the same. This is not a realistic scenario, but -- it is the basis of all the interactive models we have used up to now! It is commonly used as a basis for comparison in epidemiological modeling.
        """
    )
    return


@app.cell(hide_code=True)
def _(D00, DI1, N, plt, t00, tI1):
    plt.figure(figsize=(8,4))
    plt.vlines(30,ymin=-0.5,ymax=2*N,linestyles='dashed',colors='black')
    plt.plot(tI1,DI1,label='D')
    plt.scatter(t00[0:28],D00[0:28],label='D obs')
    plt.xlabel('Time (days)')
    plt.ylabel('Cumulative deaths')
    plt.xlim([0,100])
    plt.ylim([0.00001,1000])
    plt.legend()
    plt.title(f"1. Outbreak with no action")
    plt.gca()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ### 2) Wait for more information before taking a response
        Given the great uncertainty one may decide to try and gather more information before taking decisive action. Waiting will give experts time to analyze, synthesize and discuss existing data and there may be more data to use in analysis. Consider that we wait one week to get more information, where upon epidemiologists estimate parameters associated with the duration of infectiousness and how transmissible the pathogen is - estimating an $\mathcal{R_0}\approx5$, and an infection fatality rate of $\approx1\%$. However, during this time the epidemic has continued to develop and a large number of deaths have occurred. Upon learning this concerning news guidance is issued to reduce close-contact interactions (at the second vertical dashed line in the plot below).
        """
    )
    return


@app.cell(hide_code=True)
def _(D00, DI2, N, plt, t00, tI2):
    plt.figure(figsize=(8,4))
    plt.vlines(30,ymin=-0.5,ymax=2*N,linestyles='dashed',colors='black')
    plt.vlines(30+7,ymin=-0.5,ymax=2*N,linestyles='dashed',colors='black')
    plt.plot(tI2,DI2,label='D')
    plt.scatter(t00[0:28],D00[0:28],label='D obs')
    plt.xlabel('Time (days)')
    plt.ylabel('Cumulative deaths')
    plt.xlim([0,100])
    plt.ylim([0.00001,1000])
    plt.legend()
    plt.title(f"2. Outbreak response following delayed action")
    plt.gca()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ### 3) Large gatherings
        Large gatherings have the potential to increase transmission - as many close-contact interactions can occur. These interactions could lead to increased transmission if they take place in settings with poor ventilation, and the pathogen is airborne.

        If we let large gatherings take place we might increase the contact rate between individuals $c$. We assume that we institute three days of gathering to commemorate those known to have died from the disease. After, guidance is issued to immediately limit close contact interactions. While this short period of gatherings do not appear to obviously change the dynamics associated with fatalities, they do change the disease dynamics - we will revisit later on.
        """
    )
    return


@app.cell(hide_code=True)
def _(D00, DI3, N, plt, t00, tI3):
    plt.figure(figsize=(8,4))
    plt.vlines(30,ymin=-0.5,ymax=2*N,linestyles='dashed',colors='black')
    plt.vlines(33,ymin=-0.5,ymax=2*N,linestyles='dashed',colors='black')
    plt.plot(tI3,DI3,label='D')
    plt.scatter(t00[0:28],D00[0:28],label='D obs')
    plt.xlabel('Time (days)')
    plt.ylabel('Cumulative deaths')
    plt.xlim([0,100])
    plt.ylim([0.00001,1000])
    plt.legend()
    plt.title(f"3. Outbreak with large gathering")
    plt.gca()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ### 4) Reducing close-contact interactions
        Recognizing the potential dangers from this outbreak, perhaps guidance is immediately issued to reduce close-contact interactions to try and stop more chains of transmission from forming. We see total deaths falling from those in the previous interventions.
        """
    )
    return


@app.cell(hide_code=True)
def _(D00, DI4, N, plt, t00, tI4):
    plt.figure(figsize=(8,4))
    plt.vlines(30,ymin=-0.5,ymax=2*N,linestyles='dashed',colors='black')
    plt.plot(tI4,DI4,label='D')
    plt.scatter(t00[0:28],D00[0:28],label='D obs')
    plt.xlabel('Time (days)')
    plt.ylabel('Cumulative deaths')
    plt.xlim([0,100])
    plt.ylim([0.00001,1000])
    plt.legend()
    plt.title(f"4. Outbreak with swift response")
    plt.gca()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ### 5) Best-case scenario
        If the only tools available are able to reduce transmission - and they are so effective that they stop all new infections, we are still left with adverse outcomes. At the point of time that decisions can be made many people are infected - these infections are already baked into the system. In the below plot we set $\beta =0$ at the time of awareness, yet we see deaths in the population continue to rise (even if there are fewer total deaths).
        """
    )
    return


@app.cell(hide_code=True)
def _(D00, DI5, N, plt, t00, tI5):
    plt.figure(figsize=(8,4))
    plt.vlines(30,ymin=-0.5,ymax=2*N,linestyles='dashed',colors='black')
    plt.plot(tI5,DI5,label='D')
    plt.scatter(t00[0:28],D00[0:28],label='D obs')
    plt.xlabel('Time (days)')
    plt.ylabel('Cumulative deaths')
    plt.xlim([0,100])
    plt.ylim([0.00001,1000])
    plt.legend()
    plt.title(f"5. Outbreak with no new infections")
    plt.gca()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## Comparing these response scenarios.

        How many lives were saved from actions taken? This is a very hard question to answer, but we can attempt to answer it with mathematical modeling. When public health actions are both fast and effective - we hope to see few adverse events; crisis averted! However, success may not always be well appreciated by society at large -- if few adverse events occurred; was the outbreak as dangerous as it was considered, and were the stringent actions taken 'really' necessary?

        Given the disease parameters we can compare how simulations with no interventions compare against those with (or even against data). Such comparisons are known as **counterfactuals** - they compare against what has not come to be. Below, we compare the interventions outlined above and calculate averted deaths by comparison with the number of deaths in the scenario in which no action was taken (Intervention 1). By this metric, we see that scenarios that took immediate action were the most effective.
        """
    )
    return


@app.cell(hide_code=True)
def _(DI1, DI2, DI3, DI4, DI5, plt):
    plt.figure(figsize=(8,4))
    figsb,axsbar = plt.subplots(2,1,squeeze=False,sharex = True
                               )
    axsbar[0,0].bar([1,2,3,4,5],[DI1[-1],DI2[-1],DI3[-1],DI4[-1],DI5[-1]])
    axsbar[1,0].bar([1,2,3,4,5],[DI1[-1]-DI1[-1],DI1[-1]-DI2[-1],DI1[-1]-DI3[-1],DI1[-1]-DI4[-1],DI1[-1]-DI5[-1]])
    axsbar[1,0].set_xlabel("Intervention")
    axsbar[0,0].set_ylabel("Total deaths")
    axsbar[1,0].set_ylabel("Deaths averted")
    plt.gca()
    return axsbar, figsb


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""So far, we have only considered deaths being observable. It may be instructive to see the full suite of population dynamics associated with each intervention. Below, we show the time-varying effective reproduction number $\mathcal{R}_e$ associated with each scenario; where the horizontal dashed line represents $\mathcal{R}_e=1$, and the vertical dashed line represents the time at which the outbreak was first reported. Below that, we show the associated population dynamics for each disease state. Note: the gathering scenario (Intervention 3) has the largest peak for  infectious individuals among these scenarios.""")
    return


@app.cell(hide_code=True)
def _(
    N,
    R_effI1,
    R_effI2,
    R_effI3,
    R_effI4,
    R_effI5,
    plt,
    tI1,
    tI2,
    tI3,
    tI4,
    tI5,
    ta,
):
    plt.figure(figsize=(8,4))
    figs,axs = plt.subplots(2,3,squeeze=False,sharex=True, sharey=True)
    axs[0,0].plot(tI1,R_effI1,label='R_e')
    axs[0,0].hlines(1,xmin=-10,xmax=ta[-1]+10,linestyles='dashed',colors='black')
    axs[0,0].vlines(30,ymin=-0.5,ymax=2*N,linestyles='dashed',colors='black')
    axs[0,1].plot(tI2,R_effI2,label='R_ei')
    axs[0,1].hlines(1,xmin=-10,xmax=ta[-1]+10,linestyles='dashed',colors='black')
    axs[0,1].vlines(30,ymin=-0.5,ymax=2*N,linestyles='dashed',colors='black')
    axs[0,2].hlines(1,xmin=-10,xmax=ta[-1]+10,linestyles='dashed',colors='black')
    axs[0,2].vlines(30,ymin=-0.5,ymax=2*N,linestyles='dashed',colors='black')
    axs[0,2].plot(tI3,R_effI3,label='R_e')
    axs[1,0].plot(tI4,R_effI4,label='R_e')
    axs[1,0].hlines(1,xmin=-10,xmax=ta[-1]+10,linestyles='dashed',colors='black')
    axs[1,0].vlines(30,ymin=-0.5,ymax=2*N,linestyles='dashed',colors='black')
    axs[1,1].plot(tI5,R_effI5,label='R_e')
    axs[1,1].hlines(1,xmin=-10,xmax=ta[-1]+10,linestyles='dashed',colors='black')
    axs[1,1].vlines(30,ymin=-0.5,ymax=2*N,linestyles='dashed',colors='black')
    axs[1,2].remove()

    axs[0,0].set_xlim([0,75])
    axs[0,1].set_xlim([0,75])
    axs[0,2].set_xlim([0,75])
    axs[1,0].set_xlim([0,75])
    axs[1,1].set_xlim([0,75])
    axs[0,0].set_ylim([0,5.2])
    axs[0,1].set_ylim([0,5.2])
    axs[0,2].set_ylim([0,5.2])
    axs[1,0].set_ylim([0,5.2])
    axs[1,1].set_ylim([0,5.2])

    axs[0,0].set_title("Intervention 1")
    axs[0,1].set_title("Intervention 2")
    axs[0,2].set_title("Intervention 3")
    axs[1,0].set_title("Intervention 4")
    axs[1,1].set_title("Intervention 5")
    axs[1,1].set_xlabel("Time (days)")

    figs.suptitle(r'Different effective reproduction numbers, $\mathcal{R}_e (t)$')
    plt.gca()
    return axs, figs


@app.cell(hide_code=True)
def _(
    DI1,
    DI2,
    DI3,
    DI4,
    DI5,
    EI1,
    EI2,
    EI3,
    EI4,
    EI5,
    II1,
    II2,
    II3,
    II4,
    II5,
    RI1,
    RI2,
    RI3,
    RI4,
    RI5,
    SI1,
    SI2,
    SI3,
    SI4,
    SI5,
    plt,
    tI1,
    tI2,
    tI3,
    tI4,
    tI5,
):
    plt.figure(figsize=(8,4))
    figpopdyn,axspd = plt.subplots(2,3,squeeze=False,sharex=True, sharey=True)
    axspd[0,0].plot(tI1,SI1,label="S")
    axspd[0,0].plot(tI1,EI1,label="E")
    axspd[0,0].plot(tI1,II1,label="I")
    axspd[0,0].plot(tI1,RI1,label="R")
    axspd[0,0].plot(tI1,DI1,label="D")
    axspd[0,0].set_xlim([0,75])
    axspd[0,0].legend()

    axspd[0,1].plot(tI2,SI2,label="S")
    axspd[0,1].plot(tI2,EI2,label="E")
    axspd[0,1].plot(tI2,II2,label="I")
    axspd[0,1].plot(tI2,RI2,label="R")
    axspd[0,1].plot(tI2,DI2,label="D")
    axspd[0,1].set_xlim([0,75])

    axspd[0,2].plot(tI3,SI3,label="S")
    axspd[0,2].plot(tI3,EI3,label="E")
    axspd[0,2].plot(tI3,II3,label="I")
    axspd[0,2].plot(tI3,RI3,label="R")
    axspd[0,2].plot(tI3,DI3,label="D")
    axspd[0,2].set_xlim([0,75])

    axspd[1,0].plot(tI4,SI4,label="S")
    axspd[1,0].plot(tI4,EI4,label="E")
    axspd[1,0].plot(tI4,II4,label="I")
    axspd[1,0].plot(tI4,RI4,label="R")
    axspd[1,0].plot(tI4,DI4,label="D")
    axspd[1,0].set_xlim([0,75])

    axspd[1,1].plot(tI5,SI5,label="S")
    axspd[1,1].plot(tI5,EI5,label="E")
    axspd[1,1].plot(tI5,II5,label="I")
    axspd[1,1].plot(tI5,RI5,label="R")
    axspd[1,1].plot(tI5,DI5,label="D")
    axspd[1,1].set_xlim([0,75])

    axspd[1,2].remove()

    axspd[0,0].set_title("Intervention 1")
    axspd[0,1].set_title("Intervention 2")
    axspd[0,2].set_title("Intervention 3")
    axspd[1,0].set_title("Intervention 4")
    axspd[1,1].set_title("Intervention 5")
    axspd[1,1].set_xlabel("Time (days)")
    axspd[0,0].set_ylabel("Population")
    axspd[1,0].set_ylabel("Population")
    plt.gca()
    return axspd, figpopdyn


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        # Takeaways
        Here, we viewed interventions through a simple lens using an epidemiological SEIRD model to simulate various scenarios of (in)action in response to learning of a disease outbreak. In particular, we focused our intervening attempts via limiting the rate at which individuals were coming into close contact with one another. We found that **(a)** in outbreaks a few deaths can very quickly become many, rising at an exponential rate; **(b)** taking prompt action can be much better than awaiting additonal information (note: given uncertainties taking several actions may be better than relying on one alone); **(c)** even if an intervention that could block new infections could be applied (Intervention 5) more deaths will follow from those who are already exposed and infected.

        In particular, (note all simulation coding is at the bottom of this notebook) the epidemiological parameters used are:
        Pre-awareness contact rate $c = 10$/day, chance an infectious contact results in transmission, $p = 0.1$, such that $\beta = c.p = 1$. The average time from transmission to disease onset was $1/\nu = 2$ days; and the average duration of infectiousness $1/\gamma = 5$ days. This disease had a basic reproduction number, $\mathcal{R}_0=5$; and an infection fatality ratio (IFR) of $f=0.01$. The disease spreads with these properties until time $t_a = 30$ days after the intial exposure in all scenarios.
        Once the guidance to limit close-contact interactions is issued, we assume that the average number of contacts reduces to $c=2$/day; and changing $\beta=0.2$. Given people do not live alone, and they may still need to partake in essential activities it is unlikely to reduce to 0 (as assumed in Intervention 5).
        In Intervention 3, we assumed that c increases to $c=16$/day for the three days of large gatherings, such that $\beta=1.6$ during this time. Large gatherings may have significant impacts by driving transmissions (e.g. [Kendall et al. 2024. Science](https://doi.org/10.1126/science.adm8103)) but including such events is often not accounted for in disease modeling.
        To keep things simple in the above we focused on changes to the contact rate, but many interventions may also target the risk that a contact can turn into a transmission event.
        """
    ).callout("info")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        # Try your own intervention
        Use the sliders below to create your own scenario(s). Assess how both the effectiveness (which may depend on strength of reduction, coverage of population, and adherence by population) and the timing of intervention are important to defining success. The first vertical dashed line represents the time of awareness; and the second vertical dashed line represents the time of intervention.
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    tc = mo.ui.slider(start=0,stop=50,step=0.1,value = 25, label="time of intervention after awareness (days)")
    c0= mo.ui.slider(start=0,stop=20,step=0.1,value = 10, label="intial contact rate (per day)")
    c1= mo.ui.slider(start=0,stop=20,step=0.1,value = 10, label="changed contact rate (per day)")
    q0= mo.ui.slider(start=0,stop=0.25,step=0.01,value = 0.1, label="intial probability of a contact resulting in transmission")
    q1= mo.ui.slider(start=0,stop=0.25,step=0.01,value = 0.1, label="changed probability of a contact resulting in transmission")

    mo.md(f"""
    {c0}
    {q0}

    {tc}

    {c1}
    {q1}
    """)
    return c0, c1, q0, q1, tc


@app.cell(hide_code=True)
def _(D00, DI1, Ea, Ia, N, Ra, Sa, plt, t00, tI1, ta, tc):
    plt.figure(figsize=(8,4))
    plt.vlines(30,ymin=-0.5,ymax=2*N,linestyles='dashed',colors='black')
    plt.vlines(30+tc.value,ymin=-0.5,ymax=2*N,linestyles='dashed',colors='black')
    plt.plot(ta,Sa,label='S')
    plt.plot(ta,Ea,label='E')
    plt.plot(ta,Ia,label='I')
    plt.plot(ta,Ra,label='R')
    plt.plot(tI1,DI1,label='D')
    plt.scatter(t00[0:28],D00[0:28],label='D obs')
    plt.xlabel('Time (days)')
    plt.ylabel('Cumulative deaths')
    plt.xlim([0,100])
    plt.ylim([0.00001,max(Sa)+10])
    #plt.yscale("log")
    plt.legend()
    plt.title(f"Interactive outbreak with intervention")
    plt.gca()
    return


@app.cell(hide_code=True)
def _(D00i, Da, N, plt, t00, ta, tc):
    plt.figure(figsize=(8,4))
    plt.vlines(30,ymin=-0.5,ymax=2*N,linestyles='dashed',colors='black')
    plt.vlines(30+tc.value,ymin=-0.5,ymax=2*N,linestyles='dashed',colors='black')
    plt.plot(ta,Da,label='D')
    plt.scatter(t00[0:28],D00i[0:28],label='D obs')
    plt.xlabel('Time (days)')
    plt.ylabel('Deaths')
    plt.xlim([0,tc.value+50])
    plt.ylim([0.00001,max(Da)+10])
    plt.legend()
    plt.gca()
    return


@app.cell(hide_code=True)
def _(N, R_eff, plt, ta, tc):
    plt.figure(figsize=(8,4))
    plt.hlines(1,xmin=-10,xmax=ta[-1]+10,linestyles='dashed',colors='black')
    plt.vlines(30,ymin=-0.5,ymax=2*N,linestyles='dashed',colors='black')
    plt.vlines(30+tc.value,ymin=-0.5,ymax=2*N,linestyles='dashed',colors='black')
    plt.plot(ta,R_eff,label='R_e')
    plt.xlim([0,100])
    plt.ylim([0,max(R_eff)*1.1])
    plt.ylabel(r"$\mathcal{R}_e(t)$")
    plt.xlabel("Time (days)")
    plt.gca()
    return


@app.cell(hide_code=True)
def _(DNAi, Da, plt):
    plt.figure(figsize=(8,4))
    figsbi,axsbari = plt.subplots(2,1,squeeze=False,sharex = True
                               )
    axsbari[0,0].bar(["No action","Intervention"],[DNAi[-1],Da[-1]])
    axsbari[1,0].bar(["No action","Intervention"],[DNAi[-1]-DNAi[-1],round(DNAi[-1]-Da[-1],3)])
    axsbari[1,0].scatter(["No action","Intervention"],[DNAi[-1]-DNAi[-1],round(DNAi[-1]-Da[-1],3)],color="blue")
    axsbari[1,0].set_xlabel("Intervention")
    axsbari[0,0].set_ylabel("Total deaths")
    axsbari[1,0].set_ylabel("Deaths averted")
    plt.gca()
    return axsbari, figsbi


@app.cell
def _():
    def SEIRDmodel(u,t,beta,nu,gamma,f,N):
        S, E, I, R, D = u
        dSdt = -beta * S * I/ N
        dEdt = beta * S * I/N - nu * E
        dIdt = nu * E - gamma * I
        dRdt = (1-f) * gamma * I 
        dDdt = f * gamma * I
        return dSdt, dEdt, dIdt, dRdt, dDdt
    return (SEIRDmodel,)


@app.cell
def _():
    #fixed quantities
    N = 100000  #population size
    u0 = [N-1,1,0,0,0]
    t00 = range(31)

    nu = 1/2
    gamma = 1/5
    f= 0.01  #infection fatality ratio
    return N, f, gamma, nu, t00, u0


@app.cell
def _(
    N,
    SEIRDmodel,
    c0,
    c1,
    f,
    gamma,
    np,
    nu,
    odeint,
    q0,
    q1,
    t00,
    tc,
    u0,
):
    #parameters for outbreak intervention walk through
    beta0 = 10*0.1
    beta1 = 2*0.1
    betagather = 16*0.1

    tfirst = np.linspace(0,30,300)
    twait = np.linspace(30,35,70)
    tafter = np.linspace(35,365,1000)
    tgather = np.linspace(30,33,30)
    taftergather = np.linspace(33,365,1000)
    tclose = np.linspace(30,365,1000)

    t0 = np.linspace(0,30+tc.value,100)
    t1 = np.linspace(30+tc.value,30+tc.value+365,500)
    betainteractive0 = c0.value*q0.value
    betainteractive1 = c1.value*q1.value


    #first 30 days with sampling every day
    ret00 = odeint(SEIRDmodel, u0, t00, args=(beta0, nu, gamma, f, N))
    S00,E00,I00,R00,D00 = ret00.T
    #with more sampling for model plotting
    ret0 = odeint(SEIRDmodel, u0, t0, args=(beta0, nu, gamma, f, N))
    S0,E0,I0,R0,D0 = ret0.T

    #Intervention 1 - no action
    tI1 = np.linspace(0,365,3000)
    retInt1 = odeint(SEIRDmodel, u0, tI1, args=(beta0, nu, gamma, f, N))
    SI1,EI1,II1,RI1,DI1 = retInt1.T

    #Intervention 2 - wait then respond
    retInt2a = odeint(SEIRDmodel, u0, np.concatenate([tfirst,twait]), args=(beta0, nu, gamma, f, N))
    SI2a,EI2a,II2a,RI2a,DI2a = retInt2a.T
    retInt2b = odeint(SEIRDmodel, [SI2a[-1],EI2a[-1],II2a[-1],RI2a[-1],DI2a[-1]], tafter, args=(beta1, nu, gamma, f, N))
    SI2b,EI2b,II2b,RI2b,DI2b = retInt2b.T
    SI2 = np.concatenate([SI2a, SI2b])
    EI2 = np.concatenate([EI2a, EI2b])
    II2 = np.concatenate([II2a, II2b])
    RI2 = np.concatenate([RI2a, RI2b])
    DI2 = np.concatenate([DI2a, DI2b])
    tI2 = np.concatenate([tfirst,twait, tafter])

    #Intervention 3 - large gathering
    retInt3a = odeint(SEIRDmodel, u0, tfirst, args=(beta0, nu, gamma, f, N))
    SI3a,EI3a,II3a,RI3a,DI3a = retInt3a.T
    retInt3b = odeint(SEIRDmodel, [SI3a[-1],EI3a[-1],II3a[-1],RI3a[-1],DI3a[-1]], tgather, args=(betagather, nu, gamma, f, N))
    SI3b,EI3b,II3b,RI3b,DI3b = retInt3b.T
    retInt3c = odeint(SEIRDmodel, [SI3b[-1],EI3b[-1],II3b[-1],RI3b[-1],DI3b[-1]], taftergather, args=(beta1, nu, gamma, f, N))
    SI3c,EI3c,II3c,RI3c,DI3c = retInt3c.T
    SI3 = np.concatenate([SI3a, SI3b, SI3c])
    EI3 = np.concatenate([EI3a, EI3b, EI3c])
    II3 = np.concatenate([II3a, II3b, II3c])
    RI3 = np.concatenate([RI3a, RI3b, RI3c])
    DI3 = np.concatenate([DI3a, DI3b, DI3c])
    tI3 = np.concatenate([tfirst,tgather, taftergather])

    #Intervention 4 - immediately reduce close-contact interactions
    retInt4a = odeint(SEIRDmodel, u0, tfirst, args=(beta0, nu, gamma, f, N))
    SI4a,EI4a,II4a,RI4a,DI4a = retInt4a.T
    retInt4b = odeint(SEIRDmodel, [SI4a[-1],EI4a[-1],II4a[-1],RI4a[-1],DI4a[-1]], tclose, args=(beta1, nu, gamma, f, N))
    SI4b,EI4b,II4b,RI4b,DI4b = retInt4b.T
    SI4 = np.concatenate([SI4a, SI4b])
    EI4 = np.concatenate([EI4a, EI4b])
    II4 = np.concatenate([II4a, II4b])
    RI4 = np.concatenate([RI4a, RI4b])
    DI4 = np.concatenate([DI4a, DI4b])
    tI4 = np.concatenate([tfirst,tclose])

    #Intervention 5 - no new transmission after awareness reached
    retInt5a = odeint(SEIRDmodel, u0, tfirst, args=(beta0, nu, gamma, f, N))
    SI5a,EI5a,II5a,RI5a,DI5a = retInt5a.T
    retInt5b = odeint(SEIRDmodel, [SI5a[-1],EI5a[-1],II5a[-1],RI5a[-1],DI5a[-1]], tclose, args=(0, nu, gamma, f, N))
    SI5b,EI5b,II5b,RI5b,DI5b = retInt5b.T
    SI5 = np.concatenate([SI5a, SI5b])
    EI5 = np.concatenate([EI5a, EI5b])
    II5 = np.concatenate([II5a, II5b])
    RI5 = np.concatenate([RI5a, RI5b])
    DI5 = np.concatenate([DI5a, DI5b])
    tI5 = np.concatenate([tfirst,tclose])

    #Interactive
    #from day 0 to changepoint
    ret00interactive = odeint(SEIRDmodel, u0, t00, args=(betainteractive0, nu, gamma, f, N))
    S00i,E00i,I00i,R00i,D00i = ret00interactive.T
    ret0 = odeint(SEIRDmodel, u0, t0, args=(betainteractive0, nu, gamma, f, N))
    S0,E0,I0,R0,D0 = ret0.T
    u1 = [S0[-1],E0[-1],I0[-1],R0[-1],D0[-1]]
    ret1 = odeint(SEIRDmodel, u1, t1, args=(betainteractive1, nu, gamma, f, N))
    S1,E1,I1,R1,D1 = ret1.T
    #no action
    retNAinteractive = odeint(SEIRDmodel, u0, np.concatenate([t00,[100,365]]), args=(betainteractive0, nu, gamma, f, N))
    SNAi,ENAi,INAi,RNAi,DNAi = retNAinteractive.T


    Sa = np.concatenate([S0, S1])
    Ea = np.concatenate([E0, E1])
    Ia = np.concatenate([I0, I1])
    Ra = np.concatenate([R0, R1])
    Da = np.concatenate([D0, D1])
    ta = np.concatenate([t0, t1])

    R_eff =np.concatenate([S0/N*betainteractive0/gamma, S1/N*betainteractive1/gamma])

    R_effI1 = SI1/N*beta0/gamma
    R_effI2 = np.concatenate([SI2a/N*beta0/gamma, SI2b/N*beta1/gamma])
    R_effI3 = np.concatenate([SI3a/N*beta0/gamma, SI3b/N*betagather/gamma, SI3c/N*beta1/gamma])
    R_effI4 = np.concatenate([SI4a/N*beta0/gamma, SI4b/N*beta1/gamma])
    R_effI5 = np.concatenate([SI5a/N*beta0/gamma, SI5b/N*0/gamma])
    return (
        D0,
        D00,
        D00i,
        D1,
        DI1,
        DI2,
        DI2a,
        DI2b,
        DI3,
        DI3a,
        DI3b,
        DI3c,
        DI4,
        DI4a,
        DI4b,
        DI5,
        DI5a,
        DI5b,
        DNAi,
        Da,
        E0,
        E00,
        E00i,
        E1,
        EI1,
        EI2,
        EI2a,
        EI2b,
        EI3,
        EI3a,
        EI3b,
        EI3c,
        EI4,
        EI4a,
        EI4b,
        EI5,
        EI5a,
        EI5b,
        ENAi,
        Ea,
        I0,
        I00,
        I00i,
        I1,
        II1,
        II2,
        II2a,
        II2b,
        II3,
        II3a,
        II3b,
        II3c,
        II4,
        II4a,
        II4b,
        II5,
        II5a,
        II5b,
        INAi,
        Ia,
        R0,
        R00,
        R00i,
        R1,
        RI1,
        RI2,
        RI2a,
        RI2b,
        RI3,
        RI3a,
        RI3b,
        RI3c,
        RI4,
        RI4a,
        RI4b,
        RI5,
        RI5a,
        RI5b,
        RNAi,
        R_eff,
        R_effI1,
        R_effI2,
        R_effI3,
        R_effI4,
        R_effI5,
        Ra,
        S0,
        S00,
        S00i,
        S1,
        SI1,
        SI2,
        SI2a,
        SI2b,
        SI3,
        SI3a,
        SI3b,
        SI3c,
        SI4,
        SI4a,
        SI4b,
        SI5,
        SI5a,
        SI5b,
        SNAi,
        Sa,
        beta0,
        beta1,
        betagather,
        betainteractive0,
        betainteractive1,
        ret0,
        ret00,
        ret00interactive,
        ret1,
        retInt1,
        retInt2a,
        retInt2b,
        retInt3a,
        retInt3b,
        retInt3c,
        retInt4a,
        retInt4b,
        retInt5a,
        retInt5b,
        retNAinteractive,
        t0,
        t1,
        tI1,
        tI2,
        tI3,
        tI4,
        tI5,
        ta,
        tafter,
        taftergather,
        tclose,
        tfirst,
        tgather,
        twait,
        u1,
    )


if __name__ == "__main__":
    app.run()
