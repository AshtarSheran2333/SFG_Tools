module SFG_TYPES
	
	!reader could work with special BOXDATA string that can be parsed... (something like format string, can replace gro file in case of trr)
	!flags to analyze:
	
	!W -number- (water molecule)
	!I -number- -name- -total number of atoms- -which atom is center of mass- (ion) - needs some thinking in case of multiatomic ones...
	!OH -number- -name- -number of OHs (IDK what about aluminalike ones)
	!S -number- number of atoms to be skipped
	!FOH -number- free OH-
	
	!                                                       EXAMPLE:
	
	!                           W 100 I 2 oneatom 1 1 OH 4 silanol 1 OH 2 geminal 2 I 3 multi 6 3 S 60
	!   file contains 100 watermolecules,
	!   followed by 2 single atom ions,
	!   followed by 4 metal-o-h,
	!   followed by 2 metal-2O-h,
	!   followed by 3 multiatomic ions consiting of 6 atoms 3rd atom is taken as center of mass,
	!   skipping 60 atomic positions (can be solid)
	
	!then reader can hold array of waters, array of ion groups, array of hydroxyl groups and array of free OH-
	!later each of the groups can have different parameters to calculate A and M...
	
	!this will be a lot of work implementing
	
	!or maybe even someting like index file... this will definitely need to go in v2
	
	type atom
		real(8), dimension(3) ::                            position
		real(8), dimension(3) ::                            velocity
	end type atom
	
	type water
		type(atom) ::                                       o
		type(atom) ::                                       h1
		type(atom) ::                                       h2
	end type water

end module SFG_TYPES
