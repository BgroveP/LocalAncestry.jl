
using PkgTemplates

t = Template(;
           user="BgroveP",
           authors=["Bjarke G. Poulsen", "Emre Karaman"],
           plugins=[
               License(name="MIT"),
               Git(),
               GitHubActions(),
           ],
           dir="."
       )

t("AlleleOrigins")