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


    mo.md(
        f"""
        # Extended infectious disease dynamics models
        #### Dr. Stephen Beckett, University of Maryland (beckett@umd.edu)
        **Welcome** to the simulation exercise section for Day 4 of BSCI439C /    BIOL708F: "Infectious disease dynamics: a systems approach". Well done on opening the second simulation notebook! 🙌  

        Recall, the key point in working through this interactive notebook is not in learning how to code, which is challenging and goes beyond what we offer in this course, but to see how concepts we have discussed can be translated into models and simulations -- and the types of characteristics and dynamical outcomes these models describe.

        Do not hesitate to ask Stephen and Gabi questions, we do not want you to get stuck. 🧗🚧🙋

        There are 2 different model exercises to work through in this interactive notebook. Enjoy!

        ## About this notebook
        This notebook was designed and written using a marimo.io notebook which is coded in python - a general computing programming language. The bonus of marimo is that it can be used as an interactive notebook environment with reactivity! This means that we can use objects such as sliders to change input parameters, that will automatically update outputs - such as figures showing simulation data. This makes the notebook easy to use -- even without needing to learn computer programming! **Remember**, you can download images by right clicking on them, or finding the 'Export output as PNG' in the ... of the appropriate figure cell.
        """
    )
    return mo, np, odeint, plt


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        #1) What's up with a little delay? The SEIR model

        As discussed this week, one of the assumptions built into the differential equation models such as the SIR and SIRS models we studied in week 1 is that rates are characterized by an exponential distribution i.e., $\beta$ the transmission rate describes an average duration of flowing from $S$ to $I$ that follows an exponential distribution -- that is; mechanistically this model suggests that many individuals become near instantly infectious once they are contacted by an infected individual. Conversely, it also suggests some individuals take a very long time for the onset of infection to begin. Clearly this is a few steps away from what we might expect to see in reality for several infectious diseases e.g., respiratory viral infections where we may first expect the virus to proliferate in the host before the host becomes infectious. One method to account for the delay between transmission and infectiousness is to add one (and sometimes more) intermediate compartments. Here, we include one extra compartment $E$ to represent individuals who have been exposed to a disease, but are not yet infectious. The average time spent in the exposed compartment is given by ($1/{\nu}$). In doing so, we **(A)** change the distribution of time that individuals spend on average in the infected stages and **(B)** we change the model to distinguish between infected and infectious individuals. In the below model we will fix the average time between transmission and recovery (time spent in compartments both $E$ and $I$) as: $(1/{\nu}) + (1/{\gamma}) = 10$ days.

        The SEIR model can be written as: 
        $$\begin{align}\frac{dS}{dt} =& -\beta SI \nonumber\\
        \frac{dE}{dt} =& \beta SI - \nu E \nonumber\\
        \frac{dI}{dt} =& \nu E - \gamma I \nonumber\\
        \frac{dR}{dt} =& \gamma I\nonumber\end{align}$$

        and a coded definition is below. Here, with the below slider we can investigate how changing the relative average time spent in the $E$ and $I$ compartments affects disease dynamics.
        """
    )
    return


@app.cell
def _():
    def SEIRmodel(u,t,beta,nu,gamma):
        S, E, I, R = u
        dSdt = -beta * S * I
        dEdt = beta * S * I - nu * E
        dIdt = nu * E - gamma * I
        dRdt = gamma * I 
        return dSdt, dEdt, dIdt, dRdt
    return (SEIRmodel,)


@app.cell(hide_code=True)
def _(mo):
    nuslider = mo.ui.slider(0, 10, step = 0.1, value = 3, label='Average time (days) before onset 1/𝜈 :')

    mo.md(
        f""" 
        **Epidemiological paramters:**

        {nuslider}

        Average time (days) to recover 1/$\gamma$:$\quad$    10 - (1/𝜈)

        Transmission rate, *β*:  0.4 /day
        """
    )
    return (nuslider,)


@app.cell(hide_code=True)
def _(SEIRmodel, nuslider, odeint, plt, t):
    # Integrate the SEIR equations over the time grid, t.
    gamma = 1/(10 - nuslider.value)

    beta = 0.4
    u02 = [(1-1/10000), 1/10000, 0, 0]
    ret2 = odeint(SEIRmodel, u02, t, args=(beta, 1/nuslider.value, gamma))
    S2, E2, I2, R2 = ret2.T

    plt.figure(figsize=(8,4))
    plt.plot(t,S2,label='S')
    plt.plot(t,E2,label='E')
    plt.plot(t,I2,label='I')
    plt.plot(t,R2,label='R')
    plt.xlabel('Time (days)')
    plt.ylabel('Population fraction')
    plt.legend()
    plt.title(f"Simulation of epidemiological population dynamics with the SEIR model")
    plt.gca()
    return E2, I2, R2, S2, beta, gamma, ret2, u02


@app.cell(hide_code=True)
def _(mo):
    mo.accordion(
        {
            "**Question 1.1**: When is this simulation not well defined?  (click for answer)": mo.md(r"""This system is not well defined when 1/$\nu$ is set to either 0 or 10. Note, 1/0 is not well defined! $\newline$ At the extremes we are suggesting that either 1/$\nu$ = 0 or 1/$\gamma$ = 0 to which there is no solution."""),
            "**Question 1.2**: Do you think the mathematical definition for $\mathcal{R}_0$ is the same as that in the previously investigated SIR and SIRS models?  (click for answer)":mo.md(
        r"""
        Recall the definition for the basic reproduction number as: the average number of infections caused by one infectious individual in an otherwise susceptible population. Previously, for the SIR and SIRS models we had:
        $$\mathcal{R}_0 = \frac{\beta}{\gamma}$$
        Another way of looking at this is how many infections are generated during the time an individual is infectious.

        Here, in the SEIR model $\Large \frac{\large 1}{\large \gamma}$ still defines the average time an individual is infectious; and infectious individuals still transmit at rate $\beta$. Hence, the same definition of $\mathcal{R}_0$ applies! Later, we will look at scenarios with different basic reproduction numbers.    
        """
    ).callout("info")
                    })
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        #2) The dangers of asymptomatic transmission
        Previously, we have studied models that have had one type of transmission - all infections were considered to act in the same average way. Here, we will extend our model to consider two different types of infections - some infections are associated with symptoms and potentially severe disease that can lead to mortality. Additionally, we may expect symptomatic individuals to take actions e.g., reducing contacts, wearing masks to reduce their ability to transmit a disease onward. On the other hand, some infections may be less severe and have no noticable symptoms: asymptomatic individuals; which may still have the potential to transmit the disease onwards.
        Assembling these elements together, we will:

        - split our Infected ($I$) compartment into two compartments to consider symptomatic ($I_s$) and asymptomatic ($I_a$) infections.
        - consider that symptomatic individuals are aware of their infection and may reduce their transmission rate by a factor of $\delta$.
        - assume that new infections are symptomatic a fraction $p$ of the time, and asymptomatic a fraction $(1-p)$ of the time.
        - that a proportion $f$ of symptomatic infections lead to mortality.

        In doing so, we create a model slightly different than those we have worked previously. Here, flows into one compartment can be routed into two compartments ($E \rightarrow I_a$; and $E \rightarrow I_s$); and we also have flows from two compartments into one compartment ($I_a \rightarrow R$; and $I_s \rightarrow R$). Putting all this together into a model, we find the following set of equations:

        $$\large \begin{align} \frac{dS}{dt} =& \overset{\small transmission}{\overset{\small from\ asymptomatic}{\overbrace{-β_aSI_a}}}\qquad \overset{\small transmission}{\overset{\small from\ symptomatic}{\overbrace{-β_s(1-\delta)SI_s}}} \nonumber \\
         \frac{dE}{dt} =& \overset{\small new\ transmission}{\overbrace{S\left(\beta_aI_a + \beta_s(1-\delta)I_s\right)}}  - \overset{\small onset}{\overbrace{\nu E}}  \nonumber \\
         \frac{dI_s}{dt} =& \overset{\small symptomatic\ onset}{\overbrace{(1-p)\nu E}} - \overset{\small loss\ of\ infectiousness}{\overbrace{\gamma_s I_s}}  \nonumber \\
        \frac{dI_a}{dt} =& \overset{\small asymptomatic\ onset}{\overbrace{p\nu E}} - \overset{\small recovery\ from\ I_a}{\overbrace{\gamma_a I_a}}  \nonumber \\
         \frac{dR}{dt} =& \overset{\small recovery\ from\ I_s}{\overbrace{(1-f)\gamma_s I_s}} + \overset{\small recovery\ from\ I_a}{\overbrace{\gamma_a I_a}}  \nonumber  \\
         \frac{dD}{dt} =& \overset{\small death\ due\ to\ severe\ symptoms}{\overbrace{f\gamma_s I_s}}  \nonumber  \end{align}$$
        """
    )
    return


@app.cell
def _():
    def AsymptomaticModel(u, t, beta_a, beta_s, nu, gamma_a, gamma_s, delta, p, f):
            S, E, Ia, Is, R, D = u  #S- susceptibles, E- exposed, Ia - asymptomatic, Is - symptomatic, R - recovered, D- dead
            dSdt = -beta_a * S * Ia - beta_s * ( 1 - delta ) * S * Is
            dEdt = beta_a * S * Ia + beta_s * ( 1 - delta ) * S * Is - nu * E
            dIadt = p * nu * E - gamma_a * Ia
            dIsdt = (1 - p) * nu * E - gamma_s * Is
            dRdt = gamma_a * Ia + (1 - f) * gamma_s * Is
            dDdt = f * gamma_s * Is
            return dSdt, dEdt, dIadt, dIsdt, dRdt, dDdt
    return (AsymptomaticModel,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(
         f"""
        **Question:** can you draw out the compartments and the flows between them for this system?
         """
     ).callout("success")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        """
        Next, we will define our scenario. We will consider a population of 10,000 individuals one of whom has been exposed to a virus in an otherwise susceptible population. We will consider how an outbreak might develop over the next year. 

        The disease is characterized as follows.
        """
    )
    return


@app.cell
def _(AsymptomaticModel, delta_slide, np, odeint):
    # Total population, N.
    N = 10000
    # Initially infected and recovered individuals at time 0, I0 and R0.
    E0,Ia0,Is0,R0,D0 = 0, 1/N, 0, 0, 0
    # Initial susceptible population at time 0.
    S0 = 1 - E0 - Ia0 - Is0 - R0 - D0

    # Initial conditions vector
    u0 = S0, E0, Ia0, Is0, R0, D0

    #parameters
    beta_s = 0.8
    beta_a = 0.75*beta_s
    delta = delta_slide.value
    nu = 1/2
    gamma_a = 1/5
    gamma_s = 1/5
    f = 0.01

    #simulation timing vector
    t = np.linspace(0, 365, 200) # creates a vector from 0 to 200, with 200 elements

    # Integrate the SIR equations over the time grid, t.
    # let's look across differing levels of asymptomaticity
    retAm1 = odeint(AsymptomaticModel, u0, t, args=(beta_a, beta_s, nu, gamma_a, gamma_s, delta, 0, f))
    retAm2 = odeint(AsymptomaticModel, u0, t, args=(beta_a, beta_s, nu, gamma_a, gamma_s, delta, 0.1, f))
    retAm3 = odeint(AsymptomaticModel, u0, t, args=(beta_a, beta_s, nu, gamma_a, gamma_s, delta, 0.2, f))
    retAm4 = odeint(AsymptomaticModel, u0, t, args=(beta_a, beta_s, nu, gamma_a, gamma_s, delta, 0.3, f))
    retAm5 = odeint(AsymptomaticModel, u0, t, args=(beta_a, beta_s, nu, gamma_a, gamma_s, delta, 0.4, f))
    retAm6 = odeint(AsymptomaticModel, u0, t, args=(beta_a, beta_s, nu, gamma_a, gamma_s, delta, 0.5, f))
    retAm7 = odeint(AsymptomaticModel, u0, t, args=(beta_a, beta_s, nu, gamma_a, gamma_s, delta, 0.6, f))
    retAm8 = odeint(AsymptomaticModel, u0, t, args=(beta_a, beta_s, nu, gamma_a, gamma_s, delta, 0.7, f))
    retAm9 = odeint(AsymptomaticModel, u0, t, args=(beta_a, beta_s, nu, gamma_a, gamma_s, delta, 0.8, f))
    retAm10 = odeint(AsymptomaticModel, u0, t, args=(beta_a, beta_s, nu, gamma_a, gamma_s, delta, 0.9, f))
    retAm11 = odeint(AsymptomaticModel, u0, t, args=(beta_a, beta_s, nu, gamma_a, gamma_s, delta, 1, f))
    S, E, Ia, Is, R, D = retAm5.T #our example sim!
    return (
        D,
        D0,
        E,
        E0,
        Ia,
        Ia0,
        Is,
        Is0,
        N,
        R,
        R0,
        S,
        S0,
        beta_a,
        beta_s,
        delta,
        f,
        gamma_a,
        gamma_s,
        nu,
        retAm1,
        retAm10,
        retAm11,
        retAm2,
        retAm3,
        retAm4,
        retAm5,
        retAm6,
        retAm7,
        retAm8,
        retAm9,
        t,
        u0,
    )


@app.cell(hide_code=True)
def _(mo):
    delta_slide = mo.ui.slider(0,1,step=0.05,value = 0.8,label="Level of transmission reduction by symptomatic individuals ẟ:")
    loglin = mo.ui.dropdown(options={"linear":"linear", "log":"log"},
                            value="linear", # initial value
                            label="Choose scaling to show on plot y-axis")

    mo.md(
        f"""
        **Epidemic parameters**

        Symptomatic individuals may take some actions e.g., self-isolating, wearing masks, etc.  that reduce their ability to transmit the infection onward. Here $\delta$ = 0 suggests symptomatic individuals do not take any steps that could reduce their transmission, while $\delta$ = 1 suggests actions taken by symptomatic individuals reduce their transmission rate to 0.

        {delta_slide}   

        transmission rate from symptomatic individuals, βs: 0.8 /day

        transmission rate from symptomatic individuals, βa: 0.75*βs /day

        Average latent period before onset of infectiousness, (1/ν) : 2 days

        Average recovery time, with (1/γa) = (1/γs): 5 days

        fraction of cases that are asymptomatic, p: 0.4

        infection fatality ratio (how many deaths per infection) f: 0.01

        {loglin}
        """
    )
    return delta_slide, loglin


@app.cell(hide_code=True)
def _(D, E, Ia, Is, R, S, loglin, plt, t):
    plt.figure(figsize=(8,4))
    plt.plot(t,S,label='S')
    plt.plot(t,E,label='E')
    plt.plot(t,Ia,label='I$_a$')
    plt.plot(t,Is,label='I$_s$')
    plt.plot(t,R,label='R')
    plt.plot(t,D,label='D')
    plt.xlabel('Time (days)')
    plt.ylabel('Population fraction')
    plt.yscale(loglin.value)
    plt.ylim([10**-4,10**0.1])
    plt.legend(loc='right')
    plt.title("Simulation of epidemiological population dynamics with the SEIaIsRD model")
    plt.gca()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        f"""
        We now have the SEIaIsRD asymptomatic simulation running. This model has 6 disease states. Can you follow the dynamics? This can be quite tricky when many of the lines are close together -- go back up and see if looking at the dynamics on a logarithmic scale could be useful! How do the dynamics change as the level of transmission reduction by symptomatic individuals is altered?

        To take a deeper dive, in the below plots we evaluate how model behavior changes with uncertainty surrounding the fraction of infections that are asymptomatic. Each point in the below corresponds to the output from one SEIaIsRD simulation.

        Note, **the y-axes are dynamic** in these plots - they respond to the range in the simulation data (which can change when you alter the interactive slider). Be careful to look at the quantitative numbers as well as the qualitative shapes in your comparisons.
        """
    ).callout("success")
    return


@app.cell(hide_code=True)
def _(
    delta_slide,
    plt,
    retAm1,
    retAm10,
    retAm11,
    retAm2,
    retAm3,
    retAm4,
    retAm5,
    retAm6,
    retAm7,
    retAm8,
    retAm9,
):
    plt.figure(figsize=(8,4))
    plt.subplot(1,3,1)
    #Here want total infections occurred in sim -- quick way is to add up everyone who got to the adsorping states R and D (and must have passed through Ia or Is). That is the last [-1] and second to last [-2] elements in the return outputs at the end of the timeseries [-1].
    Infectedlist = [retAm1.T[-2][-1]+retAm1.T[-1][-1],retAm2.T[-2][-1]+retAm2.T[-1][-1],retAm3.T[-2][-1]+retAm3.T[-1][-1],retAm4.T[-2][-1]+retAm4.T[-1][-1],retAm5.T[-2][-1]+retAm5.T[-1][-1],retAm6.T[-2][-1]+retAm6.T[-1][-1],retAm7.T[-2][-1]+retAm7.T[-1][-1],retAm8.T[-2][-1]+retAm8.T[-1][-1],retAm9.T[-2][-1]+retAm9.T[-1][-1],retAm10.T[-2][-1]+retAm10.T[-1][-1],retAm11.T[-2][-1]+retAm11.T[-1][-1]]
    #plot as a line
    plt.plot([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0],[ round(elem, 4) for elem in Infectedlist])
    #add simulation points
    plt.scatter([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0],[ round(elem, 4) for elem in Infectedlist])
    plt.xlabel("Proportion asymptomatic (p)")
    plt.ylabel("Proportion cumulatively infected")
    hh=plt.subplot(1,3,2)
    plt.title(f"SEIaIsRD model predicted mortality when symptomatic individuals reduce transmission by {100*delta_slide.value}%")
    plt.axis("off")
    plt.subplot(1,3,3)
    #now the same, but just looking at all those who entered D.
    #lines
    plt.plot([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0],[retAm1.T[-1][-1],retAm2.T[-1][-1],retAm3.T[-1][-1],retAm4.T[-1][-1],retAm5.T[-1][-1],retAm6.T[-1][-1],retAm7.T[-1][-1],retAm8.T[-1][-1],retAm9.T[-1][-1],retAm10.T[-1][-1],retAm11.T[-1][-1]])
    #points
    plt.scatter([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0],[retAm1.T[-1][-1],retAm2.T[-1][-1],retAm3.T[-1][-1],retAm4.T[-1][-1],retAm5.T[-1][-1],retAm6.T[-1][-1],retAm7.T[-1][-1],retAm8.T[-1][-1],retAm9.T[-1][-1],retAm10.T[-1][-1],retAm11.T[-1][-1]])
    plt.xlabel("Proportion asymptomatic (p)")
    plt.ylabel("Proportion deceased")
    plt.gca()
    return Infectedlist, hh


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        f"""
        ### Question:
        What proportion of asymptomatic infections is most deadly, and how does it depend on the level of reduced symptomatic transmission?
        """).callout("success")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        This model comes from [Park et al. 2023. PNAS nexus](https://doi.org/10.1093/pnasnexus/pgad106). You just investigated the dynamics shown in Figure 2! Check and see if your system map corresponds to Fig 2A.

        **core message 1**: asymptomatic transmission can represent a double-edged sword -- each asymptomatic infection is mild and may not lead to severe consequences; but -- asymptomatic individuals may be more likely to transmit infection onwards as they are unaware they are infectious. At a population scale, many more infections can translate into an increase in the number of infections with more severe outcomes.

        **core message 2**: this model represents two transmission pathways - there are two ways to be infecious: either symptomatically, or asymptomatically; and each route is characterised by its own parameters. This introduces a serious complication for our previous understanding of the basic reproduction number $\mathcal{R}_0.$ When more than one transmission route is available we have to consider how an infection following each route might contribute to generating new infections. In this model the reproduction numbers associated with an asymptomatic infection are $\large \mathcal{R}_a = \beta_a / \gamma_a$; and for symptomatic infection as: $\large \mathcal{R}_s = \beta_s / \gamma_s$. Putting this together in the full model above we recover the basic reproduction number for the SEIaIsRD model as: 
        $$\large \mathcal{R}_0 = p\mathcal{R}_a + (1-p)(1-\delta)\mathcal{R}_s.$$
        """
    ).callout("info")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.accordion(
        {
            r"**Question 2.1:** using the epidemiological parameter values above, calculate the terms for $\mathcal{R}_0$. What happens when $\delta=0.25$?  (click for answer)" : mo.md(r""" When $\delta$ = 0.25; then we find that $\Large \mathcal{R}_a = 0.75\times0.8 \times5$ is equal to $\Large (1-\delta) \mathcal{R}_s = 0.75\times 0.8 \times 5 = 3$.
            This means we have:
            $${\Large \mathcal{R}_0 = 3p + 3(1-p)}$$
            This means that the expected new infections caused by asymptomatic individuals matches those by symptomatic individuals.

            Further, as $p X + (1-p) X = X$ (when p is between 0 and 1), this result is insensitive to the fraction of new infections that are asymptomatic. With these particular choices of epidemic parameters then, regardless of the choice of p, the epidemic strength is the same, = 3.

            This explains why the above plot showing the proportion of cumulative infections appears flat with respect to p when $\delta = 0.25$. Each of these theoretical epidemics has the same strength, indicating the same number of people become infected, despite the differences in relative numbers flowing through each transmission route, and the consequent differing level of associated severity.

            Note, this also indicates a critical threshold in the dynamics for $\delta$!"""),
            r"**Question 2.2:** Consider, how might $\delta$ act as a critical threshold?  (click for answer)" : mo.md(r"""
            **Note this is a long answer, and considers a lot of details.**.

            When $\delta = 0.25$ contributions to $\mathcal{R}_0$ from the asymptomatic and symptomatic transmission routes are equal (see last question). When $\delta>0.25$ then those with symptomatic infections are taking more actions to reduce their transmission potential; and the contributions from asymptomatic transmission to generating new infections become larger than the corresponding contributions from symptomatic transmission. When $\delta>0.25$ an asymptomatic infection is expected to generate more new infections; whereas when $\delta<0.25$ a symptomatic infection is expected to generate more new infections. This does not mean this is a threshold for there being more asymptomatic infections than symptomatic infections though, which is still governed by the fraction p.

            Another threshold appears at $\delta=0.75$; this is where the contributions to new infections from symptomatic individuals $\Large (1-\delta)\mathcal{R}_s = 1$. When $\delta >=0.75$ then in an outbreak composed only of symptomatic infection (when p=0) the outbreak has an $\mathcal{R}_0 <1$ and the outbreak cannot take off. We can take this logic further!

           If there is to be an outbreak in this simulation, then $\mathcal{R}_0>1$. Walking through some algebra we find:
            $${\Large \mathcal{R}_0 = \mathcal{R}_a p + (1-\delta)(1-p)\mathcal{R}_s>1}$$
            $${\Large \mathcal{R}_a p + (1-\delta)\mathcal{R}_s - (1-\delta)\mathcal{R}_s p>1}$$
            $${\Large (\mathcal{R}_a - (1-\delta)\mathcal{R}_s) p + (1-\delta)\mathcal{R}_s >1}$$
            $${\Large (\mathcal{R}_a - (1-\delta)\mathcal{R}_s) p >1 - (1-\delta)\mathcal{R}_s}$$
            $${\Large p > \frac{1 - (1-\delta)\mathcal{R}_s}{(\mathcal{R}_a - (1-\delta)\mathcal{R}_s)}}$$
            Solving with our epidemiological parameter values for $\beta$, $\gamma$ we have a condition for an outbreak, dependent on p and \delta being:
            $${\Large p > \frac{1- 4(1-\delta)}{3 - 4(1-\delta)}}$$
            What happens when $\delta =1$?
            $${\Large p > \frac{1}{3}}$$
            An outbreak will only occur when p>1/3. Check, does this match the simulation output above?        
            """),
            r"**Question 2.3:** What happens to the SEIaIsRD model when p = 0?  (click for answer)" : mo.md(r"""
            When p=1, all infections are asymptomatic and lead to recovery. With the model and parameters chosen here, this means that no infections go via the sympotmatic compartment, and there are no fatal infections. So, in this case when p=1, the SEIaIsRD model could be written as a SEIaR model, which has the same properties as the SEIR model.""")
        }
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## **Note on the basic reproduction number:**

    The SEIaIsRD model has two transmission pathways which means the basic reproduction number cannot be simply calculated as $\beta /\gamma$. If there are characteristically different types of infection with their own properties, then both the relative amount of each type of infection, as well as the expected number of transmission events within the window an infected individual is infectious need to be accounted for.

    The basic reproduction number can also be more complicated when there are multiple ways to leave infected compartments. In all our modeling work so far, we have only considered the effects of a focal disease on population dynamics. What if we consider an additional term to account for all other types of mortality that might be occuring in a population (e.g., traffic accidents, stroke, etc.)? Considering the SEIR model framework, if we assume all-cause mortality is independent of disease status and occurs at rate m, and additionally consider a population birth rate r to balance mortality, we get this set of equations:

    The SEIR model with birth and all-cause mortality can be written as: 
    $$\begin{align}\frac{dS}{dt} =& r -\beta SI - mS \nonumber\\
    \frac{dE}{dt} =& \beta SI - \nu E - mE\nonumber\\
    \frac{dI}{dt} =& \nu E - \gamma I - mI\nonumber\\
    \frac{dR}{dt} =& \gamma I - mR \nonumber\end{align}$$

    Here, some individuals who are infected (move into E) die before they become infectious (I). Additionally, some individuals who are infectious (I) die before they recover. We have to account for these factors in our calculation of $\mathcal{R}_0$ - the average number of expected new infections caused by one infected individual in an otherwise susceptible population. Here:
    $$\Large \mathcal{R}_0 = \overset{\large fraction\ going}{\overset{\large from\ E\rightarrow I}{\overbrace{(\frac{\nu}{\nu + m})}}}\quad.\quad\overset{\large average}{\overset{\large time\ infectious}{\overbrace{(\frac{1}{\gamma + m})}}}.\quad \overset{\large transmission}{\overset{\large rate}{\overbrace{\beta}}}.$$

    """).callout("info")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        f"""
        # Fin. {mo.icon("fluent-emoji-flat:1st-place-medal")}

        You completed this interactive notebook exploring some additional epidemiology models {mo.icon("fluent-emoji-flat:clapping-hands")}. In total, you have examined:

        * SIR model {mo.icon("fluent-emoji-flat:1st-place-medal")}
        * SIRS model {mo.icon("fluent-emoji-flat:1st-place-medal")}
        * stochastic SIR model {mo.icon("fluent-emoji-flat:1st-place-medal")}
        * SEIR model {mo.icon("fluent-emoji-flat:1st-place-medal")}
        * SEIaIsRD model {mo.icon("fluent-emoji-flat:1st-place-medal")}

        You should have gained an understanding of some of the principles behind designing epidemiological models, an appreciation for the richness of dynamics that even these simple epidemiology models can display; how considering the basic and effective reproduction number can help explain the observed dynamics in large populations - and can be helpful to think about how/what interventions might be useful. The utility of these simple models is that we can analyze them and isolate how different flows/interactions might change population-level disease dynamics. This becomes much harder as the system we are considering, and the number of elements and flows, grow. While we have considered several different models, these models do not reflect the full suite of modeling for all infectious disease. Different diseases, and different disease contexts/interventions, have differing characteristics which means they may need particular flows and elements that have not been showcased here. Additionally, there are additional modeling techniques (e.g. network models, agent based models) we did not discuss that could be adopted to model infectious disease spread.

        If you have extra time, consider investigating the interactive dashboard corresponding to the modeling that Dr. Mallory Harris presented on Tuesday: 

        [http://mallory-harris.shinyapps.io/divided-disease](http://mallory-harris.shinyapps.io/divided-disease)

        """
    )
    return


if __name__ == "__main__":
    app.run()
