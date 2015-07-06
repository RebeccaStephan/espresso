#include "parser.hpp"

#ifdef REACTION_ENSEMBLE
#include "reaction_ensemble.hpp"

int tclcommand_reaction_ensemble_print_status(Tcl_Interp *interp){
	char buffer[TCL_DOUBLE_SPACE];
	if(current_reaction_system.nr_single_reactions == 0){
		sprintf(buffer,"Reaction System is not initialized\n");
		Tcl_AppendResult(interp, buffer, (char *)NULL);
	}else{
		sprintf(buffer,"Reaction System is the following:\n");
		Tcl_AppendResult(interp, buffer, (char *)NULL);
		sprintf(buffer,"Volume %f\n", current_reaction_system.volume);
		Tcl_AppendResult(interp, buffer, (char *)NULL);
		for(int single_reaction_i=0;single_reaction_i<current_reaction_system.nr_single_reactions;single_reaction_i++){
			sprintf(buffer, "#Reaction %d# \n", single_reaction_i);
			Tcl_AppendResult(interp, buffer, "educt types:\n", (char *)NULL);
			//educt types and coefficients
			for(int i=0;i<current_reaction_system.reactions[single_reaction_i]->len_educt_types;i++){
				sprintf(buffer, "%d " , current_reaction_system.reactions[single_reaction_i]->educt_types[i]);
				Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
			}
			sprintf(buffer, "\neduct coefficients: \n");
			Tcl_AppendResult(interp, buffer, (char *)NULL);
			for(int i=0;i<current_reaction_system.reactions[single_reaction_i]->len_educt_types;i++){
				sprintf(buffer, "%d " , current_reaction_system.reactions[single_reaction_i]->educt_coefficients[i]);
				Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
			}
				
    			Tcl_AppendResult(interp, "\n", (char *)NULL);
			//product types and coefficients
			Tcl_AppendResult(interp, "product types:\n", (char *)NULL);
			for(int i=0;i<current_reaction_system.reactions[single_reaction_i]->len_product_types;i++){
				sprintf(buffer, "%d " , current_reaction_system.reactions[single_reaction_i]->product_types[i]);
				Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
			}
			sprintf(buffer, "\nproduct coefficients: \n");
			Tcl_AppendResult(interp, buffer, (char *)NULL);
			for(int i=0;i<current_reaction_system.reactions[single_reaction_i]->len_product_types;i++){
				sprintf(buffer, "%d " , current_reaction_system.reactions[single_reaction_i]->product_coefficients[i]);
				Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
			}				
			//reaction_constant
			sprintf(buffer, "\nequilibrium constant: %f " , current_reaction_system.reactions[single_reaction_i]->equilibrium_constant);
			Tcl_AppendResult(interp, buffer, "\n", (char *)NULL);
		}
		
	}
	return TCL_OK;
}

int tclcommand_add_reaction(Tcl_Interp *interp, int argc, char **argv){
	single_reaction* new_reaction =(single_reaction *)malloc(sizeof(single_reaction));
	
	argc -= 1; argv += 1;
	if(ARG1_IS_S("equilibrium_constant")){
		ARG_IS_D(2, new_reaction->equilibrium_constant);
		argc-=2; argv+=2;
	}
	if(ARG1_IS_S("educt_types")){
		argc-=1; argv+=1;
		int next_added_type;
		int* educt_types=(int*)malloc(sizeof(int));
		int educt_type_counter=0;
		while(ARG_IS_I(1,next_added_type)){
			argc-=1; argv+=1;
			educt_types[educt_type_counter]=next_added_type;
			educt_type_counter+=1;
			realloc(educt_types,sizeof(int)*(educt_type_counter+1));
		}
		new_reaction->len_educt_types=educt_type_counter;
		new_reaction->educt_types=educt_types;
	}else{
		return TCL_ERROR;
	}
	if(ARG1_IS_S("educt_coefficients")){
		argc-=1; argv+=1;
		int next_added_type_coeff;
		int* educt_coefficients=(int*) malloc(sizeof(int));
		int educt_type_counter=0;
		while(ARG_IS_I(1,next_added_type_coeff)){
			argc-=1; argv+=1;
			educt_coefficients[educt_type_counter]=next_added_type_coeff;
			educt_type_counter+=1;
			realloc(educt_coefficients,sizeof(int)*(educt_type_counter+1));
		}
		new_reaction->educt_coefficients=educt_coefficients;
	}else{
		return TCL_ERROR;
	}
	if(ARG1_IS_S("product_types")){
		argc-=1; argv+=1;
		int next_added_type;
		int* product_types=(int*)malloc(sizeof(int));
		int product_type_counter=0;
		while(ARG_IS_I(1,next_added_type)){
			argc-=1; argv+=1;
			product_types[product_type_counter]=next_added_type;
			product_type_counter+=1;
			realloc(product_types,sizeof(int)*(product_type_counter+1));
		}
		new_reaction->len_product_types=product_type_counter;
		new_reaction->product_types=product_types;
	}else{
		return TCL_ERROR;
	}
	if(ARG1_IS_S("product_coefficients")){
		argc-=1; argv+=1;
		int next_added_type_coeff;
		int* product_coefficients=(int*) malloc(sizeof(int));
		int product_type_counter=0;
		while(ARG_IS_I(1,next_added_type_coeff)){
			argc-=1; argv+=1;
			product_coefficients[product_type_counter]=next_added_type_coeff;
			product_type_counter+=1;
			realloc(product_coefficients,sizeof(int)*(product_type_counter+1));
			new_reaction->product_coefficients=product_coefficients;
			//check for terminus of string
			if(argc<2) {
				break;
			}
		}
	new_reaction->nu_bar=calculate_nu_bar(new_reaction->educt_coefficients, new_reaction-> len_educt_types,  new_reaction->product_coefficients, new_reaction->len_product_types);
	
	}else{
		return TCL_ERROR;
	}
	
	//if everything is fine:
	current_reaction_system.reactions[current_reaction_system.nr_single_reactions]=new_reaction;
	current_reaction_system.nr_single_reactions+=1;
	
	//assign different types an index in a growing list that starts at and is incremented by 1 for each new type
	update_type_index(new_reaction->educt_types, new_reaction->len_educt_types, new_reaction->product_types, new_reaction->len_product_types); 
	return TCL_OK;
}

int tclcommand_reaction_ensemble(ClientData data, Tcl_Interp *interp, int argc, char **argv){
	int err = TCL_OK;
	if(argc == 1){
		err= tclcommand_reaction_ensemble_print_status(interp);
	} else {
		if(current_reaction_system.reactions==NULL && ARG1_IS_S("volume")){
			// reaction_ensemble volume <volume> must be called first, it also initializes the reaction_system
			create_current_reaction_system_struct();
			double volume;
			ARG_IS_D(2,volume);
			current_reaction_system.volume=(float) volume;
		}
		if(ARG1_IS_S("do")) {
			do_reaction();
		}
		if( ARG1_IS_S("add_reaction") ) {
			//add reaction. note all calls of this function have to be done before the first call of "initialize()"
			tclcommand_add_reaction(interp, argc, argv);
		
		}
		if(ARG1_IS_S("initialize")) {
			//initializes particle types in part gc, has to be called after all reactions are added
			initialize();
		}
		
		if(ARG1_IS_S("set_default_charge_of_type")) {
			//needs to be called for each type individually after initialize was called
			int type;
			ARG_IS_I(2,type);
			double charge;
			ARG_IS_D(3,charge);
			current_reaction_system.charges_of_types[find_index_of_type(type)]=charge;
			
		}
		
		if( ARG1_IS_S("free_memory")) {
			free_reaction_ensemble();
		}
		
		if( ARG1_IS_S("set_water_type")) {
			//check for warter_type for making autodissociation of water possible in C
			ARG_IS_I(2,current_reaction_system.water_type);
		}
		
	}
		
	return gather_runtime_errors(interp,err);
}



#endif
