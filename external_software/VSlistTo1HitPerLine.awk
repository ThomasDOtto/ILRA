# VSlistTo1HitPerLine by Paul Kitts.
# Converts the vecscreen text list output to one line giving the coordinates
# for each vector segment in the format:
# VecScreen_Category   ID_string   start_position   end_position
# The default is to report Strong, Moderate, and Weak matches and also segments 
# of Suspect Origin. Reporting of any category can be suppressed by including 
# suspect=0, weak=0, moderate=0 or strong=0 on the command line. 
# "No hits" will be reported for any Query sequence that had no matches in any 
# of the selected categories, unless no_hits=0 is included on the command line.
# VecScreen errors will be reported unless errors=0 is included on the command line.
# Usage:
# VSlistTo1HitPerLine.awk [supect=0] [weak=0]... [no_hits=0] vecscreen_results_file

BEGIN {
 strong=1; moderate=1; weak=1; suspect=1;   # default categories to report
 no_hits=1; errors=1;                       # defaults is to report both
}


!/^[0-9 \t]+$/ { hits_to_report = ""}

hits_to_report {
  printf("VecScreen_%-8s\t%s\t%s\n", hits_to_report, ID, $0);
  hits++;
  next;
}

/^>Vector / {
  if (ID) {                         # check if previous Query had no hits or errors
    if (error_found) printf("VecScreen_ERROR   \t%s\n", ID);
    else if (hits == 0 && no_hits) printf("VecScreen_No_Hits \t%s\n", ID);
  }
  ID = $2;
  hits = error_found = 0;
  next;
}

strong && /^Strong/ { hits_to_report = "Strong";  next }

moderate && /^Moderate/ { hits_to_report = "Moderate";  next }

weak && /^Weak/ { hits_to_report = "Weak";  next }

suspect && /^Suspect/ { hits_to_report = "Suspect";  next }

errors && /^ERROR|WARNING/ { error_found++; next }


END {
  if (ID) {                         # check if last Query had no hits or errors
    if (error_found) printf("VecScreen_ERROR   \t%s\n", ID);
    else if (hits == 0 && no_hits) printf("VecScreen_No_Hits \t%s\n", ID);
  }
}

