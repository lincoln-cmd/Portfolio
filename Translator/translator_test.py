import googletrans

translator = googletrans.Translator()

code = googletrans.LANGCODES
code_key = code.keys()
code_value = code.values()

text = ''

while True:
    if text == '':
        print('\033[95m' + "Insert text :" + '\033[0m', end = ' ')
        text = str(input())
    
    else:
        print('\033[94m' + "Insert language or code you want to translate. If you want to search or see the list of the language, input 'language' or 'langcode' :" + '\033[0m', end = ' ')
        dst = str(input())
        
        if dst.lower() == 'language':
            print('\033[92m' + "language : " + '\033[0m', code_key)
            
        elif dst.lower() == 'langcode':
            print('\033[92m' + "language & code : " + '\033[0m', code)
            
        elif dst.lower() not in code_key and dst.lower() not in code_value:
            print('\033[91m' + "Unavailable language. Please check the correct language." + '\033[0m')
            continue
            
        else:
            if dst.lower() in code_key:
                dst = code_key[dst]
            result = translator.translate(text, dest = dst)
            
            print('\033[96m' + f"{result.text}" + '\033[0m')
            print('\033[96m' + f"{result.pronunciation}" + '\033[0m')
            text = ''
            