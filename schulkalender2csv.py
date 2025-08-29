# Dieses Skript extrahiert die Termine aus den (in Hamburg üblichen) Excel-Schulkalendern,
# in denen jeweils 6 Monate nebeneinander angeordnet sind und darunter in 30 Zeilen die
# einzelnen Tage mit den Terminen.
# Jede Monatsspalte umfasst dabei 4 Excel-Spalten: 1. Tag, 2. Wochentag, 3. Termine, 4. KW (an Montagen)
# Mehrere Termine an einem Tag sind durch Semikolon getrennt.
# Darunter folgt das gleiche für die Monate Februar-Juli.
#
# Die Ausgabe enthät 2 Spalten im CSV Format:
# 1. Datum (yyyy-mm-dd), 2. Terminbeschreibung

import os, openpyxl

def extract_texts_with_positions(file_path):
    """
    Extracts all text values from an Excel file (.xlsx) 
    and returns a list of dictionaries with sheet, cell, row, column, and value.
    """
    wb = openpyxl.load_workbook(file_path, data_only=True)
    results = []

    for sheet_name in wb.sheetnames:
        ws = wb[sheet_name]
        for row in ws.iter_rows():
            for cell in row:
                if cell.value not in (None, ""): # and cell.value.strip():
                    value = cell.value
                    # Wenn der Wert ein float ist, der eine ganze Zahl darstellt, in int umwandeln
                    if isinstance(value, float) and value.is_integer():
                        value = int(value)
                    results.append({
                        "sheet": sheet_name,
                        "cell": cell.coordinate,   # e.g. "B12"
                        "row": cell.row,
                        "col": cell.column,       # integer index
                        "value": str(value).replace("\n"," ")
                    })
    return results

if __name__ == "__main__":
    file_name = "" # hier ggf. Dateinamen eingeben
    if file_name == "":
        file_name = input("xlsx-Datei: (z.B. Downloads\schulkalender.xlsx): ")
    file_path = os.path.expandvars("%USERPROFILE%\\")+file_name
    texts = extract_texts_with_positions(file_path)

    months = ["Januar", "Februar", "März", "April", "Mai", "Juni", "Juli",
              "August", "September", "Oktober", "November", "Dezember"]
    start_col=[1,5,9,13,17,21]
    month_col={}
    skip_rows=[]
    events=[]
    for entry in texts:
        row,col,value = entry['row'],entry['col'],entry['value']
        if value.split(" ")[0] == "Schulkalender":
            skip_rows.append(row)
        if row in skip_rows:
            continue
        if col-1 in start_col:
            weekday = value
            continue
        if (col in start_col and value.split(" ")[0] in months):
            m,y=value.split(" ")
            mnr = months.index(m)+1
            if mnr<10: mnr="0"+str(mnr)
            month_col[col]=f"{y}-{mnr}"
            continue
        if col in start_col:
            day=int(value)
            if day<10: day="0"+str(day)
            continue
        if col-2 in start_col:
          values = value.split(";")
          for v in values:
            events.append([f"{month_col[col-2]}-{day}", v.strip()])

    events.sort(key = lambda e: e[0])

    for e in events:
        print('"'+e[0]+'","'+e[1].replace('"','""')+'"')
