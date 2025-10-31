class TeamAsparagine:
    def __init__(self, members):
        self.members = members

    def show_members(self):
        print("Team Members:")
        for member in self.members.keys():
            print(f"- {member}")
        print()

    def show_usernames(self):
        print("Slack Usernames:")
        for member, info in self.members.items():
            print(f"{member}: {info['Slack Username']}")
        print()

    def show_countries(self):
        print("Countries:")
        for member, info in self.members.items():
            print(f"{member}: {info['Country']}")
        print()

    def show_hobbies(self):
        print("Hobbies:")
        for member, info in self.members.items():
            print(f"{member}: {info['Hobby']}")
        print()

    def show_affiliations(self):
        print("Affiliations:")
        for member, info in self.members.items():
            print(f"{member}: {info['Affiliation']}")
        print()

    def show_fav_sequences(self):
        print("Favourite DNA Sequences:")
        for member, info in self.members.items():
            print(f"{member}: {info['DNA Sequence of Favourite Gene']}")
        print()



member_team_asparagine = {
    "Ayşenur Akcan": {
        "Slack Username": "@Ayşenur Akcan",
        "Country": "Türkiye",
        "Hobby": "Mending something",
        "Affiliation": "Istanbul Technical University",
        "DNA Sequence of Favourite Gene": "TP3- ATGGAGGAGCCGCAGTCAGATCCTAGCGTCGAGCCCCCTCTGAGTCAGGAAACATTTTCAGACCTATGGAAACTACTTCCTGAAAACAACGTTCT"
    },
    "Oni David Sunday": {
        "Slack Username": "@Sunday",
        "Country": "Nigeria",
        "Hobby": "Reading",
        "Affiliation": "Wageningen University and Research",
        "DNA Sequence of Favourite Gene": "AKT1 - TCTTCTCCCACTCCCCAGCAAAGTCCCCCTTTTGTGAGTGTAGCTGCCAGTACCTAGGTGAATGGTTGACTCCCCTCGGAGCCTTCCT..."
    },
    "Noor Ul Ain Amir": {
        "Slack Username": "@Noor Ul Ain Amir",
        "Country": "Pakistan",
        "Hobby": "Baking",
        "Affiliation": "Forman Christian College, University",
        "DNA Sequence of Favourite Gene": "SHH - ACAAGCTCTCCAGGCTTGCTACCATTTAAAATCAGACTCTTTTTGTCTTTTGATTGCTGTCTCGCGACCCAACTCCGATGTGTTCCGTTATCAGCGGCCGGCAGCCTGCCATTCCAGCCCCT..."
    },
    "Ozlem Kalkan": {
        "Slack Username": "@anyavala...",
        "Country": "Ukraine- Turkey",
        "Hobby": "doing sport",
        "Affiliation": "Bonn University",
        "DNA Sequence of Favourite Gene": "8OHY - ATGGCTCAAAGCACTGTTTTGCCTATGCACTGTTTATACGGGATATTCCTCGAAGGCAACCTAAAGATTCAAAAGAATGATCAGGAAGGACTAAAGAAGTTTAAAGATAATATCAAAAAATT…"
    },
    "Daria Kriuchkova": {
        "Slack Username": "@Daria",
        "Country": "Ukraine",
        "Hobby": "Knitting",
        "Affiliation": "Kharkiv Polytechnic Institute",
        "DNA Sequence of Favourite Gene": "BRCA1 - ATGGAAGTTGTCATTTTGTGTTTCCAGGATTTATTTGCTCTTCGTGTCTTTGGGTAGCTGG…"
    },
    "Minenhle Mayisela": {
        "Slack Username": "@Minenhle Mayisela",
        "Country": "Eswatini and South Africa",
        "Hobby": "Reading",
        "Affiliation": "University of the Witwatersrand",
        "DNA Sequence of Favourite Gene": "LMNA - GACAAATTCCTTGACCCGAGGAGGATAGGGATGTGGCCTTCGGTCTTTCCTCGCAGCTCCGGGGCAAGCTAGGAGTGGGATGGAAGTCGAGGTCCCTAATTTTTTAAGGGGAGGGTGCGGGGAGAAGGGGTAGTATGCGGAAACAGAGCGGGTATGAAGCTGGCTAACGCCGCGCGCCCCCTCCCAGGACCCGCTCCTGCCCCGCGCCGG"
    },
    "Seventh Member": {
        "Slack Username": "@...",
        "Country": "...",
        "Hobby": "...",
        "Affiliation": "...",
        "DNA Sequence of Favourite Gene": "..."
    },
    "Eighth Member": {
        "Slack Username": "@...",
        "Country": "...",
        "Hobby": "...",
        "Affiliation": "...",
        "DNA Sequence of Favourite Gene": "..."
    }
}



if __name__ == "__main__":
    team = TeamAsparagine(member_team_asparagine)
    team.show_members()
    team.show_usernames()
    team.show_countries()
    team.show_hobbies()
    team.show_affiliations()
    team.show_fav_sequences()
