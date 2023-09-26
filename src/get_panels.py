from panelapp import Panelapp, api, queries
import requests, json

panel_output = open("panelapp_response.json", "w")

api_url = "https://panelapp.genomicsengland.co.uk/api/docs/?format=openapi"
response = requests.get(api_url)
json.dump(response.json(), panel_output, indent=2)

panel_output.close()


# panels = queries.get_all_panels()

# for panel, panel_info in panels.items():
#     print(panel, panel_info, file=panel_output)
#     break

# panel_output.close()
