import pandas as pd

def calculate_residence_times(interactions_per_frame, frame_time):
    # Now stores {lipid_id: (start_time, protein_resid)}
    active_interactions = {} 
    all_events = []          
    
    print("Calculating residence times and binding sites...")
    
    for step_index, frame_data in enumerate(interactions_per_frame):
        # current_lipids is now a dict: {lipid_id: protein_resid}
        current_lipids = frame_data['lipids'] 
        current_time = step_index * frame_time
        
        # 1. Check for ended interactions
        for lipid in list(active_interactions.keys()):
            if lipid not in current_lipids:
                start_time, prot_resid = active_interactions.pop(lipid)
                duration = current_time - start_time
                
                all_events.append({
                    'lipid_id': lipid,
                    'protein_resid': prot_resid, # <--- NEW!
                    'start_time': start_time,
                    'end_time': current_time,
                    'duration': duration
                })
                
        # 2. Check for new interactions
        for lipid, prot_resid in current_lipids.items():
            if lipid not in active_interactions:
                # Record the time and the binding site
                active_interactions[lipid] = (current_time, prot_resid)
                
    # 3. Handle lipids still interacting at the end
    final_time = (len(interactions_per_frame) - 1) * frame_time
    for lipid, (start_time, prot_resid) in active_interactions.items():
        all_events.append({
            'lipid_id': lipid,
            'protein_resid': prot_resid, # <--- NEW!
            'start_time': start_time,
            'end_time': final_time,
            'duration': final_time - start_time
        })
        
    df = pd.DataFrame(all_events)
    if not df.empty:
        df = df.sort_values(by=['lipid_id', 'start_time']).reset_index(drop=True)
        
    return df
    