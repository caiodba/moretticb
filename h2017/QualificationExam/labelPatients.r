## LABEL LIST
.labels=c(
	"esquerdo",
	"esquerdo",
	"esquerdo", 
	"direito",
	"esquerdo",
	"direito",
	"direito",
	"direito",
	"esquerdo", 
	"esquerdo",
	"esquerdo",
	"direito",
	"direito",
	"direito",
	"esquerdo",
	"esquerdo",
	"direito",
	"esquerdo",
	"direito",
	"direito",
	"direito",
	"direito",
	"esquerdo",
	"esquerdo",
	"direito",
	"direito", 
	"esquerdo",
	"esquerdo",
	"esquerdo",
	"esquerdo",
	"direito", 
	"direito",
	"direito",
	"direito",
	"esquerdo"
)

## ALIAS LIST
.aliases=c(
	"P1",
	"P2",
	"P3",
	"P4",
	"P5",
	"P6",
	"P7",
	"P8",
	"P9",
	"P10",
	"P11",
	"P12",
	"P12",
	"P13",
	"P14",
	"P15",
	"P16",
	"P17",
	"P13",
	"P29",
	"P7",
	"P30",
	"P31",
	"P5",
	"P32",
	"P16",
	"P33",
	"P34",
	"P35",
	"P11",
	"P36",
	"P37",
	"P38",
	"P12",
	"P35"
)

getLabelFromID <- function(id){
	type = substr(id,start=1,stop=1)
	if(type == "P"){ #patient
		return(.labels[match(id,.aliases)])
	} else { #healthy individual
		return("error")
	}
}
