# import web app modules
import streamlit as st
from numpy import random as rd
from plotly.graph_objs import *
from plotly.subplots import make_subplots

st.set_page_config(layout="wide")


def allele_simulation(initA, fAA, fAa, faa, pop, gen,sim):
# simulating the data    
    freqA_list_sim = [] # creating list to store A allele frequencies. Will be a list of lists
    
    for _ in range(sim):
        freqA_list = [initA] 
    
        for _ in range(gen):
            last_gen_freqA = freqA_list[-1] #last generation's frequency of A
            last_gen_freqa = 1 - last_gen_freqA #last generation's frequency of a
        
            M = fAA*last_gen_freqA**2 + fAa*2*last_gen_freqA*last_gen_freqa + faa*last_gen_freqa**2 # used to standardize the frequencies later

            #calculate the frequency of A in the current generation, accounting for finite population size
            AA = rd.binomial(pop, (fAA*last_gen_freqA**2)/M)
            Aa = rd.binomial(pop, (fAa * 2*last_gen_freqA*last_gen_freqa)/M)
            aa = rd.binomial(pop, (faa*last_gen_freqa**2)/M)
            tot = AA+Aa+aa
            freqA = (AA + 0.5*Aa)/tot
    
            freqA_list.append(freqA) # appending current frequency of A to freqA_list
    
        freqA_list_sim.append(freqA_list) # appending entire list with every generation of A frequency to the overall freqA_list_sim

# making the output plots

    fig = make_subplots(rows=1, cols=2)
    simulations = []
    
    for index, value in enumerate(freqA_list_sim):
        legend_title = f"Simulation {index + 1}"
        simulation = Scatter(x=list(range(len(value))), y=value, mode='lines', name=legend_title)
        simulations.append(simulation)
        
    for simulation in simulations:
        fig.add_trace(simulation, row=1, col=2)
        
    lastgen_freqA = []
    
    for sim in freqA_list_sim:
        lastgen_freqA.append(sim[-1])
        
    fig.add_trace(Histogram(x=lastgen_freqA, nbinsx= 10, xbins=dict(start=0, end=1.1, size = .1), showlegend=False, hoverinfo='skip', marker_line_width=1,marker_line_color="black"), row=1, col=1)
        
    fig.update_yaxes(title_text="Allele A Frequency", range=[0, 1], showline = True, linecolor = 'black', linewidth = 1, row=1, col=2)
    fig.update_xaxes(title_text="Generation", showline = True, linecolor = 'black', linewidth = 1, row=1, col=2)
    fig.update_yaxes(title_text="Count", showline = True, linecolor = 'black', linewidth = 1, row=1, col=1)
    fig.update_xaxes(title_text="Allele A Frequency", range=[0, 1.1], showline = True, linecolor = 'black', linewidth = 1, row=1, col=1)
    
    fig.update_layout(plot_bgcolor = 'white',
                  title_text="Simulation Outputs",
                  height=450, 
                  width=1000)
    return fig

left_column, right_column = st.columns([1,2])

with left_column:
    sim = st.slider('Number of Simulations:',min_value=1,max_value=100,step=1)
    fAA = st.slider('Fitness of AA:',min_value=0.,max_value=1.,step=0.05,value=1.)
    fAa = st.slider('Fitness of Aa:',min_value=0.,max_value=1.,step=0.05,value=1.)
    faa = st.slider('Fitness of aa:',min_value=0.,max_value=1.,step=0.05,value=1.)
    pop = st.select_slider('Population Size:',[10,50,100,500,1000])
    gen = st.slider("Number of Generations:",min_value=100,max_value=1000,step=100)
    initA = st.slider("Starting frequency of A:", min_value=0.01,max_value=0.99,step=0.01, value=0.5)

with right_column:
    st.header("Pop Gen Simulation & Histogram")
    graph = allele_simulation(initA, fAA, fAa, faa, pop, gen,sim)
    st.plotly_chart(graph)


