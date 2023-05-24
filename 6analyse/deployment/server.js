const express = require('express')
const basicAuth = require('basic-auth-connect')
const fs = require("fs")
const sharp = require("sharp")
const app = express()

const USERNAME = process.env.BASIC_USER_NAME
const PASSWORD = process.env.BASIC_PASSWORD

app.use(basicAuth(function (user, password) {
    const bool = USERNAME == user && PASSWORD == password
    return bool
}))

app.get("/resized/:kind/:file", (req, res) => {
    const filepath = `public/graphs/${req.params.kind}/${req.params.file}`
    console.log(filepath)
    fs.readFile(filepath, (err, data) => {
        if (err != null) {
            res.sendStatus(404)
            return
        }
        width = req.query.width ? parseInt(req.query.width, 10) : 200
        height = req.query.height ? parseInt(req.query.height, 10) : 200
        sharp(data).resize({ width, height }).toBuffer().then(resizedData => {
            res.set("Content-Type", "image/png")
            res.send(resizedData)
        })
    })
})

app.use(express.static('public'))

const port = process.env.PORT || 3000
app.listen(port, () => {
    console.log(`App listening on port ${port}.`)
})
