import plotly.express as px
import seaborn as sns

def lipid_protein_scatter_plot(df):
    fig = px.scatter(df, x="Residue ID", y="Frame", color="Lipid Type",
                     hover_data=["Protein", "Residue Name", "Lipid ID"])
    fig.show()

def lipid_protein_line_plot(df):
    df1 = df.groupby(['Frame', 'Lipid Type']).count()['Protein'].reset_index()
    fig = px.line(df1, x='Frame', y='Protein', color='Lipid Type',
                  labels={'Protein': 'Number of Interactions', 'Frame': 'Time (Frame)'})
    fig.show()

# def lipid_protein_scatter_plot(df):
#     fig = px.scatter(df, x='Residue ID', y='Lipid ID', color='Lipid Type',
#                      hover_data=['Protein', 'Residue Name', 'Frame'],
#                      labels={'Lipid ID': 'Lipid Interactions', 'Residue ID': 'Residue'})
#     fig.show()

def lipid_protein_bar_plot(df, frame=1):
    df1 = df[df['Frame'] == frame]
    df1 = df1.groupby(['Frame', 'Lipid Type']).count()['Protein'].reset_index()
    fig = px.bar(df1, x='Lipid Type', y='Protein', color='Lipid Type', labels={'Protein': 'Number of Interactions'})
    fig.show()

# def lipid_protein_histogram(df):
#     fig = px.histogram(df, x='Lipid ID', color='Lipid Type',
#                        nbins=30, marginal='box',
#                        hover_data=['Protein', 'Residue Name', 'Residue ID', 'Frame'],
#                        labels={'Lipid ID': 'Number of Interactions'})
#     fig.show()

def lipid_protein_box_plot(df, first_f=1, last_f=1):
    df1 = df[df['Frame'] > first_f]
    df1 = df1[df1['Frame'] < last_f]
    df1 = df1.groupby(['Frame', 'Lipid Type']).count()['Protein'].reset_index()
    fig = px.box(df1, x='Lipid Type', y='Protein', labels={'Protein': 'Number of Interactions'})
    fig.show()

def lipid_protein_density_heatmap(df):
    df1 = df.groupby(['Residue ID', 'Lipid ID']).count()['Frame'].reset_index()
    fig = px.density_heatmap(df1, x='Residue ID', y='Lipid ID', z='Frame', nbinsx=df1['Residue ID'].max(), nbinsy=df1['Lipid ID'].max())
    fig.show()

def lipid_protein_box_metrics(df):
    fig = px.box(df, x='Residue ID', y='Sum of all contacts', color='Lipid Type')
    fig.show()
