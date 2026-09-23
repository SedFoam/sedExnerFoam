"""
Function to read dataset from Lyn (1988)
"""


def readExpParameters(pathFileExp, caseName=None):
    caseList = []
    with open(pathFileExp, "r") as fData:
        for line in fData:
            if line[0] == "#":
                continue
            case = {}
            name, Hwater, Umean, uf, ws, dS = line[:-1].split(";")
            case["name"] = name
            case["Hwater"] = float(Hwater)
            case["Umean"] = float(Umean)
            case["uf"] = float(uf)
            case["ws"] = float(ws)
            case["dS"] = float(dS)
            if caseName == case["name"]:
                return case
            caseList.append(case)
    return caseList
