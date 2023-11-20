import os

def create_user_file():
    # Ask the user for username and password
    username = input("Enter your username: ")
    password = input("Enter your password: ")

    # Create a text file and write the username and password to it
    with open("user_XNAT.txt", "w") as file:
        file.write(f"Username: {username}\n")
        file.write(f"Password: {password}\n")

def read_user_file():
    # Read the username and password from the text file
    with open("user_XNAT.txt", "r") as file:
        lines = file.readlines()
        username = lines[0].split(":")[1].strip()
        password = lines[1].split(":")[1].strip()

    return username, password
def main():
    # Check if the file exists
    if os.path.exists("user_XNAT.txt"):
        # If the file exists, read username and password
        existing_username, existing_password = read_user_file()
        print(f"Existing Username: {existing_username}")
        print(f"Existing Password: {existing_password}")
    else:
        # If the file does not exist, create a new file and ask for username and password
        create_user_file()
        print("User file created successfully.")
        existing_username, existing_password = read_user_file()
        print(f"Existing Username: {existing_username}")
        print(f"Existing Password: {existing_password}")
    return existing_username, existing_password