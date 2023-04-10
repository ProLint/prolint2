import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

class PL2Plotter:
    def __init__(self, contacts='contacts.csv', metrics='contacts_metrics.csv'):
        self.contacts_df = pd.read_csv(contacts)
        self.metrics_df = pd.read_csv(metrics)

    def interaction_heatmap(self, metric='Sum of all contacts'):
        """
        Plots a heatmap of the interaction matrix between protein residues and lipid residues.
        """
        fig, ax = plt.subplots(figsize=(10,10)) 
        if metric not in [
            "Sum of all contacts",
            "Occupancy",
            "Longest Duration",
            "Mean Duration",
        ]:
            raise ValueError("The metric is not valid.")
        else:
            sns.heatmap(self.metrics_df.pivot_table(index='Lipid ID', columns='Residue ID', values=metric), cbar_kws={'label': metric})
            plt.title("Protein-Lipid Interaction Matrix")
            plt.show()

    # def sequence_heatmap(self, metric='Sum of all contacts'):
    #     """
    #     Plots a heatmap on the sequence of the protein.
    #     """
    #     fig, ax = plt.subplots(figsize=(10,10)) 
    #     if metric not in [
    #         "Sum of all contacts",
    #         "Occupancy",
    #         "Longest Duration",
    #         "Mean Duration",
    #     ]:
    #         raise ValueError("The metric is not valid.")
    #     else:
    #         sns.heatmap(self.metrics_df.pivot_table(index='Lipid ID', columns='Residue ID', values=metric), cbar_kws={'label': metric})
    #         plt.title("Protein-Lipid Interaction Matrix")
    #         plt.show()
        
    # def plot_occupancy(self):
    #     """
    #     Plots a bar chart showing the occupancy of each lipid type.
    #     """
    #     sns.barplot(x="Lipid Type", y="Occupancy", data=self.summary_df)
    #     plt.title("Lipid Occupancy")
    #     plt.show()
        
    # def plot_duration(self):
    #     """
    #     Plots a bar chart showing the mean duration of interactions between protein and lipid for each lipid type.
    #     """
    #     sns.barplot(x="Lipid Type", y="Mean Duration", data=self.summary_df)
    #     plt.title("Mean Interaction Duration")
    #     plt.show()
