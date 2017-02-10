#!/usr/bin/env python2

import connexion
from flask_cors import CORS

if __name__ == '__main__':
    app = connexion.App(__name__, specification_dir='./swagger/')
    app.add_api('swagger.yaml', arguments={'title': 'This REST Api allows you to submit chemical identifiers in one format and translate it into another format (e.g. SMILES -&gt; InChi)'})
    CORS(app.app)
    app.run(port=8080)
