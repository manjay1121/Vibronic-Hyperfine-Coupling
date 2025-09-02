import pandas as pd

def read_first_block_to_dataframe(filename):
    with open(filename, 'r') as file:
        # Read the first line to get N (number of data rows)
        first_line = file.readline()
        N = int(first_line.strip().split()[0])
        
        # Skip the second line (header or null bit)
        file.readline()
        
        # Read the next N lines containing 7 columns each
        data = []
        for _ in range(N):
            line = file.readline().strip()
            parts = line.split()
            # Convert each column to float (or keep as string if you prefer)
            row = [parts[0]] + [float(x) for x in parts[1:]]
#            row = [float(x) for x in parts]
            data.append(row)

    # Create pandas DataFrame with 7 columns
    column_names = [f'Col{i+1}' for i in range(7)]
    df = pd.DataFrame(data, columns=column_names)

    # Print the DataFrame
#    print(df.iloc[:, :4].to_string(index=False, header=False))

    # Return the DataFrame for follow-up calculations
    return df

def create_transformed_dataframe(df,F):
    # Create a new DataFrame with the specified columns
    df_plus = pd.DataFrame({
        'Index': df['Col1'],                       # Copy the first column as is
        'Col2_x_Col5': df['Col2'] + (df['Col5'] * F),  # Multiply 2nd and 5th columns
        'Col3_x_Col6': df['Col3'] + (df['Col6'] * F),  # Multiply 3rd and 6th columns
        'Col4_x_Col7': df['Col4'] + (df['Col7'] * F)  # Multiply 4th and 7th columns
    })
    df_minus = pd.DataFrame({
        'Index': df['Col1'],                       # Copy the first column as is
        'Col2_x_Col5': df['Col2'] - (df['Col5'] * F),  # Multiply 2nd and 5th columns
        'Col3_x_Col6': df['Col3'] - (df['Col6'] * F),  # Multiply 3rd and 6th columns
        'Col4_x_Col7': df['Col4'] - (df['Col7'] * F)  # Multiply 4th and 7th columns
    })

    return df_plus, df_minus

# Replace 'your_file.txt' with the actual filename you want to read
filename = 'coord'
F = 0.05

# Call the function and store the result in a variable
df = read_first_block_to_dataframe(filename)
df_plus, df_minus = create_transformed_dataframe(df, F)

#print("Positive displacement:")
#print(df_plus.to_string(index=False, header=False))

#print("\nNegative displacement:")
#print(df_minus.to_string(index=False, header=False))

with open("positive.xyz", "w") as f:
    f.write(df_plus.to_string(index=False, header=False))

# Save negative displacement to a file
with open("negative.xyz", "w") as f:
    f.write(df_minus.to_string(index=False, header=False))

