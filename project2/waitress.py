

from waitress import serve

serve(wsgiapp, host='0.0.0.0', port=8080)

